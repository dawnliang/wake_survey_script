function [Di, Dp, U, V, W] = calculate_drag(probeData)
    %{
        @filename   calculate_drag.m
        @author     DL + NP
        @date       July 22, 2015
        @updated    August 1, 2015 -DL
        UW 21__ NextGen wake survey dr

        calculates the parasitic & induced drag values based on the set of
        data points in a wake survey; locations of data points can be
        arbitrary.

        ** TODO: update column offsets in READ FILE section to match with final
        data formatting **
        ** TODO: recomment such that none of the test code is executed **

        based on:
            "Quantitative Three-Dimensional Low-Speed Wake Surveys" by G.W. Brune
            http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19930018258.pdf

        includes files:
            - wakeSurveyDrag.m      main executable script
            - calc2Integral.m       custom function for calculating surface
                                    integrals (over the wake)
            - calcInterpCoeff.m     custom function for interpolating data
                                    points (over the wake)

        uses data channels:
            - QC, PRESSTS, TEMPTS
            - PTP1/2, UTP1/2, VTP1/2, WTP1/2
            - TRAV_Y, TRAV_Z (TRAV_X is assumed constant)

        important equation variables
            - Y, Z: cartesian coordinates of location in wake
            - U, V, W: matrices of wake pressures
                U <=> x; tunnel axis
                V, W <=> y, z; measuring plane crossflow notation
            - ro: air density
            - p_t: total pressure
            - u_inf: freestream u-axis air velocity

        notation:
            - _inf = undisturbed freestream value
    %}

    %% READ DATA
    Cp = probeData(:,1); % unitless coefficient
    ALPHA = probeData(:,2); % degrees upflow
    BETA = probeData(:,3); % degrees crossflow
    Y = probeData(:,5); % in
    Z = probeData(:,6); % in
    QA = probeData(:,7); % psf
    PRESSTS = probeData(:,8); % psi
    TEMPTS = probeData(:,9); % deg F

    %% CALCULATE CONSTANTS & PERFORM UNIT CONVERSIONS

    %{
            - p_diff = p_t_inf - p_t (psf)
                - p_t: total pressure (local) (psi)
                - p_t_inf: freestream total pressure (psf)
            - ro: air density (slugs / ft^3)
            - u_inf: freestream u-axis air pressure (psf)
    %}

    RANKINE_CONST = 549.67;         % deg R = deg F + RANKINE CONST
    UNIVERSAL_GAS_CONSTANT = 1716;  % ft * lb / (slug * R)
    SQFT_TO_SQIN = 144;             % conversion constant

    %{
        p = ro*R*T => ro = p / (R*T)
            - p: pressure in lb/ft^2 (from PRESSTS)
            - ro: density in slugs/ft^3
            - R: universal gas constant = 1716 ft*lb/(slug*R)
            - T: temperature in Rankine (F + 459.67) (from TEMPTS)
    %}
    ro = (PRESSTS * SQFT_TO_SQIN) ./ (UNIVERSAL_GAS_CONSTANT .* (TEMPTS + RANKINE_CONST));

    %{
        p_diff = p_t_inf - p_t = Qa (1- Cp
            - p_t_inf; freestream total pressure (Qa) (psf)
            - p_t; local total pressure (Cp * Qa) (psf)
    %}
    p_diff = QA .* (1 - Cp);
    
    for i = 1:length(p_diff)
        if p_diff(i) < 0
            p_diff(i) = 0;
        end
    end

    %{
        q = 1/2 * ro*v^2 => v = sqrt(2*q/ro) = u_inf
            - q: dynamic pressure (Qa)
            - ro: density
            - v: velocity of wind
    %}
    u_inf = sqrt(2 .* QA ./ ro);

    %{
        total pressure data is given in pressure coefficients; multiply by Q &
        convert to velocity
    %}
    P = sqrt(2 .* Cp .* QA ./ ro);

    %{
        pressure data given in pressure/upflow/xflow; convert to uvw
            U = p * cos(alpha) * sin(beta)
            V = p * cos(alpha) * cos(beta)
            W = p * sin(alpha)
                where alpha = upflow, beta = xflow
    %}
    U = P .* cosd(ALPHA) .* cosd(BETA);
    V = -P .* cosd(ALPHA) .* sind(BETA);
    W = P .* sind(ALPHA);

    %% CALCULATE DRAG

    %{
        profile/parasitic drag
        D_p = double integral over the wake of 
                [(p_t_inf - p_t) + ro/2 (U* - U)(U* + U - 2U_inf)] ds

        utilises custom function calc2Integral; see documentation in function

        intermediates:
            artificial axial velocity:  U* ^2 = U^2 + 2/ro (p_t_inf - p_t)
                                        U* = sqrt(U^2 + 2/ro (p_t_inf - p_t))
    %}

    art_axial_v = sqrt(U.^2 + (2 ./ ro .* p_diff));
    U_a = (art_axial_v - U) .* (art_axial_v + U - 2 .* u_inf);
    func = p_diff + ro/2 .* U_a;

    Dp = calculate_surface_integral(Y./12, Z./12, func);

    %{
        induced drag:
        D_i = ro/2 * double integral over the wake of
                (V^2 + W^2 - u'^2) ds
            = double integral over the wake of
                ro/2 (V^2 + W^2 - u'^2) ds

        utilises custom function calc2Integral; see documentation in function

        intermediates:
            perturbation velocity (pert_v) u' = U* - U_inf
                - U* = artificial axial velocity (see profile/parasitic drag)
    %}

    pert_v = art_axial_v - u_inf;
    func = (ro ./ 2) .* (V.^2 + W.^2 - pert_v.^2);

    Di = calculate_surface_integral(Y./12, Z./12, func);
end
