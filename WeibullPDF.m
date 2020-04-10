# Vectorial function to calculate probabilities 
# from Weibull probability density function

function p = WeibullPDF(v,A,k)
  p = (k/A)*((v./A).^(k-1.0)).*exp(-(v./A).^k);
endfunction