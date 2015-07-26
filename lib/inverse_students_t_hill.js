'use strict';

// MODULES //

var erfcinv = require( 'compute-erfcinv/lib/number.js' );


// FUNCTIONS //

var pow = Math.pow,
	sqrt = Math.sqrt;

var expm1 = Math.expm1 || function(x) {
	return Math.exp(x) - 1;
};

// CONSTANTS //

var PI = Math.PI,
	ROOT_TWO = 1.414213562373095048801688724209698078e+00;

//
// The main method used is due to Hill:
//
// G. W. Hill, Algorithm 396, Student's t-Quantiles,
// Communications of the ACM, 13(10): 619-620, Oct., 1970.
//
function inverse_students_t_hill( ndf, u, pol ) {

	var a, b, c, d, q, x, y;

	if ( ndf > 1e20 ) {
		return -erfcinv( 2 * u ) * ROOT_TWO;
	}
	a = 1 / (ndf - 0.5);
	b = 48 / (a * a);
	c = ((20700 * a / b - 98) * a - 16) * a + 96.36;
	d = ((94.5 / (b + c) - 3) / b + 1) * sqrt( a * PI / 2 ) * ndf;
	y = pow(d * 2 * u, 2 / ndf);

	if (y > (0.05 + a)) {
		//
		// Asymptotic inverse expansion about normal:
		//
		x = -erfcinv( 2 * u ) * ROOT_TWO;
		y = x * x;

		if (ndf < 5) {
			c += 0.3 * ( ndf - 4.5 ) * ( x + 0.6 );
		}
		c += (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b;
		y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c - y - 3) / b + 1) * x;
		y = expm1(a * y * y, pol);
	} else {
		y = ((1 / ( ( (ndf + 6) / (ndf * y) - 0.089 * d - 0.822 ) *
		(ndf + 2) * 3 ) + 0.5 / (ndf + 4)) * y - 1) *
		(ndf + 1) / (ndf + 2) + 1 / y;
	}
	q = sqrt( ndf * y );
	return -q;
}

module.exports = inverse_students_t_hill;
