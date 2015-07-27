'use strict';

var ibeta_inv_imp = require( './ibeta_inv_imp.js' );

function betaincinv_l( p, a, b ) {
	// The argument a to the incomplete beta function inverse must be greater than zero.
	if ( a <= 0 ) {
		return NaN;
	}
	// The argument b to the incomplete beta function inverse must be greater than zero.
	if ( b <= 0 ) {
		return NaN;
	}
	// p must be inside the range [0,1]
	if( (p < 0) || (p > 1) ) {
		return NaN;
	}
	return ibeta_inv_imp( a, b, p, 1 - p )[ 0 ];
}

function betaincinv_u( p, a, b ) {
	// The argument a to the incomplete beta function inverse must be greater than zero.
	if ( a <= 0 ) {
		return NaN;
	}
	// The argument b to the incomplete beta function inverse must be greater than zero.
	if ( b <= 0 ) {
		return NaN;
	}
	// p must be inside the range [0,1]
	if( (p < 0) || (p > 1) ) {
		return NaN;
	}
	return ibeta_inv_imp( a, b, p, 1 - p )[ 1 ];
}

// EXPORTS //

module.exports = {
	'lower': betaincinv_l,
	'upper': betaincinv_u
};
