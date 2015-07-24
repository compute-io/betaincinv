'use strict';

// MODULES //

var gammaincinv = require( 'compute-gammaincinv' ),
	temme_root_finder = require( './root_finder.js'),
	newton_raphson_iterate = require( './newton_raphson.js' );


// FUNCTIONS //

var log = Math.log,
	sqrt = Math.sqrt;


//
// See:
// "Asymptotic Inversion of the Incomplete Beta Function"
// N.M. Temme
// Journal of Computation and Applied Mathematics 41 (1992) 145-157.
// Section 4.
//
function temme_method_3_ibeta_inverse( a, b, p, q ) {
	var mu,
		w, w_2, w_3, w_4, w_5, w_6, w_7, w_8, w_9, w_10,
		d, d_2, d_3, d_4,
		e1, e2, e3,
		w1, w1_2, w1_3, w1_4,
		eta0, eta,
		u, cross, lower, upper, x;
	//
	// Begin by getting an initial approximation for the quantity
	// eta from the dominant part of the incomplete beta:
	//
	if( p < q ) {
	  eta0 = gammaincinv( p, b, {
			'tail': 'upper'
	  });
 	} else {
		eta0 = gammaincinv(q, b, {
			'tail': 'lower'
		});
	}
	eta0 /= a;
	//
	// Define the variables and powers we'll need later on:
	//
	mu = b / a;
	w = sqrt(1 + mu);
	w_2 = w * w;
	w_3 = w_2 * w;
	w_4 = w_2 * w_2;
	w_5 = w_3 * w_2;
	w_6 = w_3 * w_3;
	w_7 = w_4 * w_3;
	w_8 = w_4 * w_4;
	w_9 = w_5 * w_4;
	w_10 = w_5 * w_5;
	d = eta0 - mu;
	d_2 = d * d;
	d_3 = d_2 * d;
	d_4 = d_2 * d_2;
	w1 = w + 1;
	w1_2 = w1 * w1;
	w1_3 = w1 * w1_2;
	w1_4 = w1_2 * w1_2;
	//
	// Now we need to compute the purturbation error terms that
	// convert eta0 to eta, these are all polynomials of polynomials.
	// Probably these should be re-written to use tabulated data
	// (see examples above), but it's less of a win in this case as we
	// need to calculate the individual powers for the denominator terms
	// anyway, so we might as well use them for the numerator-polynomials
	// as well....
	//
	// Refer to p154-p155 for the details of these expansions:
	//
	e1 = (w + 2) * (w - 1) / (3 * w);
	e1 += (w_3 + 9 * w_2 + 21 * w + 5) * d / (36 * w_2 * w1);
	e1 -= (w_4 - 13 * w_3 + 69 * w_2 + 167 * w + 46) * d_2 / (1620 * w1_2 * w_3);
	e1 -= (7 * w_5 + 21 * w_4 + 70 * w_3 + 26 * w_2 - 93 * w - 31) * d_3 / (6480 * w1_3 * w_4);
	e1 -= (75 * w_6 + 202 * w_5 + 188 * w_4 - 888 * w_3 - 1345 * w_2 + 118 * w + 138) * d_4 / (272160 * w1_4 * w_5);

	e2 = (28 * w_4 + 131 * w_3 + 402 * w_2 + 581 * w + 208) * (w - 1) / (1620 * w1 * w_3);
	e2 -= (35 * w_6 - 154 * w_5 - 623 * w_4 - 1636 * w_3 - 3983 * w_2 - 3514 * w - 925) * d / (12960 * w1_2 * w_4);
	e2 -= (2132 * w_7 + 7915 * w_6 + 16821 * w_5 + 35066 * w_4 + 87490 * w_3 + 141183 * w_2 + 95993 * w + 21640) * d_2  / (816480 * w_5 * w1_3);
	e2 -= (11053 * w_8 + 53308 * w_7 + 117010 * w_6 + 163924 * w_5 + 116188 * w_4 - 258428 * w_3 - 677042 * w_2 - 481940 * w - 105497) * d_3 / (14696640 * w1_4 * w_6);

	e3 = -((3592 * w_7 + 8375 * w_6 - 1323 * w_5 - 29198 * w_4 - 89578 * w_3 - 154413 * w_2 - 116063 * w - 29632) * (w - 1)) / (816480 * w_5 * w1_2);
	e3 -= (442043 * w_9 + 2054169 * w_8 + 3803094 * w_7 + 3470754 * w_6 + 2141568 * w_5 - 2393568 * w_4 - 19904934 * w_3 - 34714674 * w_2 - 23128299 * w - 5253353) * d / (146966400 * w_6 * w1_3);
	e3 -= (116932 * w_10 + 819281 * w_9 + 2378172 * w_8 + 4341330 * w_7 + 6806004 * w_6 + 10622748 * w_5 + 18739500 * w_4 + 30651894 * w_3 + 30869976 * w_2 + 15431867 * w + 2919016) * d_2 / (146966400 * w1_4 * w_7);
	//
	// Combine eta0 and the error terms to compute eta (Second eqaution p155):
	//
	eta = eta0 + e1 / a + e2 / (a * a) + e3 / (a * a * a);
	//
	// Now we need to solve Eq 4.2 to obtain x.  For any given value of
	// eta there are two solutions to this equation, and since the distribtion
	// may be very skewed, these are not related by x ~ 1-x we used when
	// implementing section 3 above.  However we know that:
	//
	//  cross < x <= 1       ; iff eta < mu
	//          x == cross   ; iff eta == mu
	//     0 <= x < cross    ; iff eta > mu
	//
	// Where cross == 1 / (1 + mu)
	// Many thanks to Prof Temme for clarifying this point.
	//
	// Therefore we'll just jump straight into Newton iterations
	// to solve Eq 4.2 using these bounds, and simple bisection
	// as the first guess, in practice this converges pretty quickly
	// and we only need a few digits correct anyway:
	//
	if ( eta <= 0 ) {
		eta = Number.MIN_VALUE;
	}
	u = eta - mu * log(eta) + (1 + mu) * log(1 + mu) - mu;
	cross = 1 / (1 + mu);
	lower = eta < mu ? cross : 0;
	upper = eta < mu ? 1 : cross;
	x = (lower + upper) / 2;
	x = newton_raphson_iterate(
		temme_root_finder(u, mu), x, lower, upper, 100
	);
	return x;
}

module.exports = temme_method_3_ibeta_inverse;
