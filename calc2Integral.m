%{
    FUNCTION R = CALC2INTEGRAL(X, Y, data)

    @filename   calc2Integral.m
    @author     DL + NP
    @date       July 22, 2015
    @updated    August 1, 2015 -DL

    calculates the double integral over a set of data. uses delaunay
    triangulation to create sets of triangles to integrate over, then
    calculates the integral using interpolation over normalised coordinates

    INPUT
        X, Y: the X & Y coordinates of the set of data points as column
            vectors
                X = [x1;    Y = [y1;	point 1: (x1,y1)
                     x2;         y2;    point 2: (x2,y2)
                     ...         ...
                     xm]         xm]    point m: (xm,ym)

        data: the values of the data to be integrated corresponding to the
            xy locations given by X & Y
                data = [d1;     data at point 1 (x1,y1)
                        d2;     data at point 2 (x2,y2)
                        ...
                        dm]     data at point m (xm,ym)

    OUTPUT
        r = the value of the double integral of the data over the set of
            points (X,Y)

    IMPLEMENTATION DETAILS

    delaunay triangulation: generates 3-point sets of delaunay triangles.
        the circumcircle of each triangle contains no other points in the
        data set.

        returns a matrix, with each row corresponding to the vertices of a
        triangle. there is no specific order to the triangles.

            delaunay(x,y) = [t11 t12 t13;   triangle 1
                             t21 t22 t23;   triangle 2
                             t31 t32 t33;   triangle 3
                             ... ... ...
                             tm1 tm2 tm3;   triangle m

        documentation:
            http://www.mathworks.com/help/matlab/ref/delaunay.html?searchHighlight=delaunay

    double integral: analytically evaluated using interpolation for the value of Q at
        any point on the triangle, which is normalised using change of
        variables (by NP in Mathematica)

        change of variables:
                  o
                / |                        m
        y     /   |     === === >          |
        |   o-----|                 (0, 1) |o
        |                                  || \
        | ___ ___ __x                      |o_-_o_____l
                                    (0, 0)      (1, 0)
        
        interpolation:
        Q = sum from i = 1:3 of Ni(x,y)*Qi
          = N1(x,y)*Q1 + N2(x,y)*Q2 + N3(x,y)*Q3

        double integral on triangle S of sum from i = 1:3 of Ni(x,y)*Qi
        = 1/6 * (Q1 + Q2 + Q3 * ( b3(c1 - 3 a2 c1 - c2 + 3 a1 c2)
            + 3 a3 (b2 c1 - b1 c2)
            + c3 (-b1 + 3 a2 b1 + b2 - 3 a1 b2))
            / (b2 c1 - b1 c2))
%}

function r = calc2Integral(X, Y, data)
    % // TEST CODE - to test the function, uncomment the code below & replace
    % //             X/Y with x/y
    
    %{
    x = gallery('uniformdata',[3600,1],1);
    y = gallery('uniformdata',[3600,1],2);
    data = gallery('uniformdata',[3600,1],3);
    %}
    
    % triangulation
    triangles = delaunay(X, Y);
    triplot(triangles, X, Y); % plot triangles
    
    r = 0;
    for i = 1:size(triangles, 1)
        % isolate the triangle
        tri = [X(triangles(i,:)) Y(triangles(i,:))];
        q = data(triangles(i,:));
        
        % interpolation & CoV
        coeff = calcInterpCoeff(tri);
        % coeff = [a1 a2 a3;
        %          b1 b2 b3;
        %          c1 c2 c3]
        
        % evaluation of double integral
        % R = 1/6 * (Q1 + Q2 + Q3 * ( b3(c1 - 3 a2 c1 - c2 + 3 a1 c2) + 3 a3 (b2 c1 - b1 c2) + c3 (-b1 + 3 a2 b1 + b2 - 3 a1 b2)) / (b2 c1 - b1 c2))
        r = r + 1/6 * (q(1) + q(2) +...
            q(3) .* (coeff(2,3)*(coeff(3,1) - 3*coeff(1,2)*coeff(3,1) - coeff(3,2) + 3*coeff(1,1)*coeff(3,2))...
            + 3*coeff(1,3)*(coeff(2,2)*coeff(3,1) - coeff(2,1)*coeff(3,2))...
            + coeff(3,3)*(-coeff(2,1)+3*coeff(1,2)*coeff(2,1)+coeff(2,2)-3*coeff(1,1)*coeff(2,2)))...
            /(coeff(2,2)*coeff(3,1) - coeff(2,1)*coeff(3,2)));
    end
end
