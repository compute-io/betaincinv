'use strict';

// MODULES //

var betainc = require( 'compute-betainc/lib/numbers.js' );

// FUNCTIONS //

var abs = Math.abs;

function ibeta_roots( a, b, target, invert ) {
	invert = invert || false;

	return function roots( x ) {
		var f1, f, f2, f1p = {},
			y;
		y = 1 - x;
		f = betainc( a, b, x, undefined, invert, true, f1p ) - target;
		f1 = f1p.value;
		if ( invert ) {
			f1 = -f1;
		}
		if (y === 0) {
			y = Number.MIN_VALUE * 64;
		}
		if (x === 0) {
			x = Number.MIN_VALUE * 64;
		}

		f2 = f1 * (-y * a + (b - 2) * x + 1);
		if( abs( f2 ) < y * x * Number.MAX_VALUE ) {
			f2 /= (y * x);
		}
		if ( invert ) {
			f2 = -f2;
		}
		// make sure we don't have a zero derivative:
		if( f1 === 0 ) {
			f1 = (invert ? -1 : 1) * Number.MIN_VALUE * 64;
		}
		return [ f, f1, f2 ];
	};
}

module.exports = ibeta_roots;
