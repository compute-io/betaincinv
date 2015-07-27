'use strict';

// MODULES //

var evaluate_polynomial = require( './evaluate_polynomial.js' ),
	zeros = require( 'compute-zeros' ),
	erfcinv = require( 'compute-erfcinv' ),
	temme_root_finder = require( './root_finder.js'),
	newton_raphson_iterate = require( './newton_raphson.js' );


// FUNCTIONS //

var abs = Math.abs,
	exp = Math.exp,
	log = Math.log,
	sqrt = Math.sqrt,
	sin = Math.sin,
	cos = Math.cos;


//
// See:
// "Asymptotic Inversion of the Incomplete Beta Function"
// N.M. Temme
// Journal of Computation and Applied Mathematics 41 (1992) 145-157.
// Section 3.
//

function temme_method_2_ibeta_inverse( z, r, theta ) {

	var eta0, s, c, terms, workspace,
		sc, sc_2, sc_3, sc_4, sc_5, sc_6, sc_7,
		co1, co2, co3, co4, co5, co6, co7, co8, co9, co10, co11, co12, co13,
		alpha, eta, s_2, c_2, u, lu,
		x;

	//
	// Get first estimate for eta, see Eq 3.9 and 3.10,
	// but note there is a typo in Eq 3.10:
	//
	eta0 = erfcinv( 2 * z );
	eta0 /= -sqrt(r / 2);

	s = sin( theta );
	c = cos( theta );
	//
	// Now we need to purturb eta0 to get eta, which we do by
	// evaluating the polynomial in 1/r at the bottom of page 151,
	// to do this we first need the error terms e1, e2 e3
	// which we'll fill into the array "terms".  Since these
	// terms are themselves polynomials, we'll need another
	// array "workspace" to calculate those...
	//
	terms = zeros( 4 );
	terms[0] = eta0;
	workspace = zeros( 6 );
	//
	// some powers of sin(theta)cos(theta) that we'll need later:
	//
	sc = s * c;
	sc_2 = sc * sc;
	sc_3 = sc_2 * sc;
	sc_4 = sc_2 * sc_2;
	sc_5 = sc_2 * sc_3;
	sc_6 = sc_3 * sc_3;
	sc_7 = sc_4 * sc_3;
	//
	// Calculate e1 and put it in terms[1], see the middle of page 151:
	//
	workspace[0] = (2 * s * s - 1) / (3 * s * c);
	co1 = [ -1, -5, 5 ];
	workspace[1] = -evaluate_polynomial( co1, s*s, 3 ) / (36 * sc_2);
	co2 = [ 1, 21, -69, 46 ];
	workspace[2] = evaluate_polynomial( co2, s*s, 4 ) / (1620 * sc_3);
	co3 = [ 7, -2, 33, -62, 31 ];
	workspace[3] = evaluate_polynomial( co3, s*s, 5 ) / (6480 * sc_4);
	co4 = [ 25, -52, -17, 88, -115, 46 ];
	workspace[4] = evaluate_polynomial( co4, s*s, 6 ) / (90720 * sc_5);
	terms[1] = evaluate_polynomial( workspace, eta0, 5 );
	//
	// Now evaluate e2 and put it in terms[2]:
	//
	co5 = [ 7, 12, -78, 52 ];
	workspace[0] = -evaluate_polynomial( co5, s*s, 4 ) / (405 * sc_3);
	co6 = [ -7, 2, 183, -370, 185 ];
	workspace[1] = evaluate_polynomial( co6, s*s, 5 ) / (2592 * sc_4);
	co7 = [ -533, 776, -1835, 10240, -13525, 5410 ];
	workspace[2] = -evaluate_polynomial( co7, s*s, 6 ) / (204120 * sc_5);
	co8 = [ -1579, 3747, -3372, -15821, 45588, -45213, 15071 ];
	workspace[3] = -evaluate_polynomial( co8, s*s, 7 ) / (2099520 * sc_6);
	terms[2] = evaluate_polynomial( workspace, eta0, 4 );
	//
	// And e3, and put it in terms[3]:
	//
	co9 = [ 449, -1259, -769, 6686, -9260, 3704 ];
	workspace[0] = evaluate_polynomial( co9, s*s, 6 ) / (102060 * sc_5);
	co10 = [ 63149, -151557, 140052, -727469, 2239932, -2251437, 750479 ];
	workspace[1] = -evaluate_polynomial( co10, s*s, 7 ) / (20995200 * sc_6);
	co11 = [ 29233, -78755, 105222, 146879, -1602610, 3195183, -2554139, 729754 ];
	workspace[2] = evaluate_polynomial( co11, s*s, 8 ) / (36741600 * sc_7);
	terms[3] = evaluate_polynomial( workspace, eta0, 3 );
	//
	// Bring the correction terms together to evaluate eta,
	// this is the last equation on page 151:
	//
	eta = evaluate_polynomial( terms, 1/r, 4 );
	//
	// Now that we have eta we need to back solve for x,
	// we seek the value of x that gives eta in Eq 3.2.
	// The two methods used are described in section 5.
	//
	// Begin by defining a few variables we'll need later:
	//
	s_2 = s * s;
	c_2 = c * c;
	alpha = c / s;
	alpha *= alpha;
	lu = (-(eta * eta) / (2 * s_2) + log(s_2) + c_2 * log(c_2) / s_2);
	//
	// Temme doesn't specify what value to switch on here,
	// but this seems to work pretty well:
	//
	if ( abs(eta) < 0.7 ) {
		//
		// Small eta use the expansion Temme gives in the second equation
		// of section 5, it's a polynomial in eta:
		//
		workspace[0] = s * s;
		workspace[1] = s * c;
		workspace[2] = (1 - 2 * workspace[0]) / 3;
		co12 = [ 1, -13, 13 ];
		workspace[3] = evaluate_polynomial( co12, workspace[0], 3 ) / (36 * s * c);
		co13 = [ 1, 21, -69, 46 ];
		workspace[4] = evaluate_polynomial( co13, workspace[0], 4 ) / (270 * workspace[0] * c * c);
		x = evaluate_polynomial( workspace, eta, 5 );
	} else {
		//
		// If eta is large we need to solve Eq 3.2 more directly,
		// begin by getting an initial approximation for x from
		// the last equation on page 155, this is a polynomial in u:
		//
		u = exp(lu);
		workspace[0] = u;
		workspace[1] = alpha;
		workspace[2] = 0;
		workspace[3] = 3 * alpha * (3 * alpha + 1) / 6;
		workspace[4] = 4 * alpha * (4 * alpha + 1) * (4 * alpha + 2) / 24;
		workspace[5] = 5 * alpha * (5 * alpha + 1) * (5 * alpha + 2) * (5 * alpha + 3) / 120;
		x = evaluate_polynomial( workspace, u, 6 );
		//
		// At this point we may or may not have the right answer, Eq-3.2 has
		// two solutions for x for any given eta, however the mapping in 3.2
		// is 1:1 with the sign of eta and x-sin^2(theta) being the same.
		// So we can check if we have the right root of 3.2, and if not
		// switch x for 1-x.  This transformation is motivated by the fact
		// that the distribution is *almost* symetric so 1-x will be in the right
		// ball park for the solution:
		//
		if ( (x - s_2) * eta < 0 ) {
			x = 1 - x;
		}
	}
	//
	// The final step is a few Newton-Raphson iterations to
	// clean up our approximation for x, this is pretty cheap
	// in general, and very cheap compared to an incomplete beta
	// evaluation.  The limits set on x come from the observation
	// that the sign of eta and x-sin^2(theta) are the same.
	//
	var lower, upper;
	if ( eta < 0 ) {
		lower = 0;
		upper = s_2;
	} else {
		lower = s_2;
		upper = 1;
	}
	//
	// If our initial approximation is out of bounds then bisect:
	//
	if( (x < lower) || (x > upper) ) {
		x = ( lower + upper ) / 2;
	}
	//
	// And iterate:
	//
	x = newton_raphson_iterate(
		temme_root_finder( -lu, alpha ), x, lower, upper, 100
	);

	return x;
}

module.exports = temme_method_2_ibeta_inverse;
