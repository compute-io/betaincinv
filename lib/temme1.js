'use strict';

// MODULES //

var polyval = require( 'compute-polynomial' ),
	zeros = require( 'compute-zeros' ),
	erfcinv = require( 'compute-erfcinv' );


// FUNCTIONS //

var exp = Math.exp,
	sqrt = Math.sqrt;


//
// See:
// "Asymptotic Inversion of the Incomplete Beta Function"
// N.M. Temme
// Journal of Computation and Applied Mathematics 41 (1992) 145-157.
// Section 2.
//

function temme_method_1_ibeta_inverse( a, b, z ) {
	var r2, eta0, terms, workspace,
		B, B_2, B_3,
		eta, eta_2,
		c, x;

	r2 = sqrt( 2 );
	//
	// get the first approximation for eta from the inverse
	// error function (Eq: 2.9 and 2.10).
	//
	eta0 = erfcinv( 2 * z );
	eta0 /= -sqrt(a / 2);

	terms = zeros( 4 );
	workspace = zeros( 7 );

	terms[ 0 ] = eta0;
	//
	// calculate powers:
	//
	B = b - a;
	B_2 = B * B;
	B_3 = B_2 * B;

	//
	// Calculate correction terms:
	//
	// See eq following 2.15:
	workspace[0] = -B * r2 / 2;
	workspace[1] = (1 - 2 * B) / 8;
	workspace[2] = -(B * r2 / 48);
	workspace[3] = -1 / 192;
	workspace[4] = -B * r2 / 3840;
	terms[1] = polyval( workspace, eta0 );
	// Eq Following 2.17:
	workspace[0] = B * r2 * (3 * B - 2) / 12;
	workspace[1] = (20 * B_2 - 12 * B + 1) / 128;
	workspace[2] = B * r2 * (20 * B - 1) / 960;
	workspace[3] = (16 * B_2 + 30 * B - 15) / 4608;
	workspace[4] = B * r2 * (21 * B + 32) / 53760;
	workspace[5] = (-32 * B_2 + 63) / 368640;
	workspace[6] = -B * r2 * (120 * B + 17) / 25804480;
	terms[2] = polyval( workspace, eta0 );
	// Eq Following 2.17:
	workspace[0] = B * r2 * (-75 * B_2 + 80 * B - 16) / 480;
	workspace[1] = (-1080 * B_3 + 868 * B_2 - 90 * B - 45) / 9216;
	workspace[2] = B * r2 * (-1190 * B_2 + 84 * B + 373) / 53760;
	workspace[3] = (-2240 * B_3 - 2508 * B_2 + 2100 * B - 165) / 368640;
	terms[3] = polyval( workspace, eta0 );
	//
	// Bring them together to get a final estimate for eta:
	//
	eta = polyval( terms, 1/a );
	//
	// now we need to convert eta to x, by solving the appropriate
	// quadratic equation:
	//
	eta_2 = eta * eta;
	c = -exp( -eta_2 / 2 );
	if ( eta_2 === 0 ) {
	  x = 0.5;
	} else {
	  x = ( 1 + eta * sqrt( (1 + c) / eta_2) ) / 2;
	}
	return x;
}

module.exports = temme_method_1_ibeta_inverse;
