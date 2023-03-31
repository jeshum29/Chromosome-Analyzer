function prob = calcSurfGoodProb(x, params)

prob = params(1) + 1/2*erf((x - params(2)) * params(3));

end