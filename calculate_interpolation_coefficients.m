%{
    FUNCTION R = CALCULATE_INTERPOLATION_COEFFICIENTS(TRI)
    
    @filename   calculate_interpolation_coefficients.m
    @author     DL + NP
    @date       July 22, 2015
    @updated    August 1, 2015 -DL

    calculates coefficients for interpolation of data at the vertices over
    a triangle

    INPUT
        tri: the (x,y) coordinates defining the vertices of the triangle to
            be interpolated over
                tri = [x1 y1;   vertex 1
                       x2 y2;   vertex 2
                       x3 y3]   vertex 3

    OUTPUT
        r: a matrix of the coefficients for interpolation with Ni(x,y);
            the ith column corresponds to the coefficients for the Ni
            interpolation
                r = [a1 a2 a3;
                     b1 b2 b3;
                     c1 c2 c3]

    DETAILS

    given:      N_i(x,y) = ai + bi*x + ci*y
                N_i(x,y) = 1
                N_j(x,y) = 0

    we can calculate ai, bi, ci using the linear systems below:
        N1(x,y):
            [ 1  = [ 1 x1 y1  [ a1
              0      1 x2 y2    b1
              0 ]    1 x3 y3 ]  c1 ]

        N2(x,y):
            [ 0  = [ 1 x1 y1  [ a2
              1      1 x2 y2    b2
              0 ]    1 x3 y3 ]  c2 ]

        N3(x,y):
            [ 0  = [ 1 x1 y1  [ a3
              0      1 x2 y2    b3
              1 ]    1 x3 y3 ]  c3 ]
%}

function r = calculate_interpolation_coefficients(tri)
     A = [1 0 0;
         0 1 0;
         0 0 1];
     locs = [ones(3,1) tri];
     r = [locs\A(:,1) locs\A(:,2) locs\A(:,3)];
end
