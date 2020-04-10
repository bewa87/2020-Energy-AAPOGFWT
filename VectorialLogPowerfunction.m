### Piecewise Defined Vectorial Power Output Function
### with Independent Wind Speed Data Variable For Five
### Parameters Logistic Regression

function pow_log = VectorialLogPowerfunction(v,z)
  n       = length(v);
  pow_log = zeros(n,1);
  for j = 1:n
      if     (v(j) >= 0 && v(j) < 2.75)
        pow_log(j) = 0.0;
      elseif (v(j) >= 2.75 && v(j) <= 12.50)
        pow_log(j) = z(1)/(z(2)+z(3)*exp(z(4)*v(j)+z(5)));
      elseif (v(j) > 12.50 && v(j) <= 25.00)
        pow_log(j) = 3110;
      elseif (v(j) > 25.00)
        pow_log(j) = 0.0;
      endif
  endfor
  pow_log = pow_log';
endfunction