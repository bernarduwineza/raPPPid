function basisValue = bspline_basis(i, degree, knots, x)
    % 

    % Recursive Cox-de Boor formula to calculate B-spline basis functions
    if degree == 0
        basisValue = (x >= knots(i)) & (x < knots(i + 1));
    else
        alpha1 = (x - knots(i)) / (knots(i + degree) - knots(i));
        alpha1(isinf(alpha1)) = 0; 
        alpha2 = (knots(i + degree + 1) - x) / (knots(i + degree + 1) - knots(i + 1));
        alpha2(isinf(alpha2)) = 0; 
        
            
        basisValue = alpha1 .* bspline_basis(i, degree - 1, knots, x) + alpha2 .* bspline_basis(i + 1, degree - 1, knots, x);
    end
end