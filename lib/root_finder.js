'use strict';

// FUNCTIONS //

var log = Math.log;

function temme_root_finder( t, a ) {
	return function roots( x ) {
		var y, f, f1, big;

		y = 1 - x;
		if (y === 0 ) {
			big = Number.MAX_VALUE / 4;
			return [ -big, -big ];
		}
		if( x === 0 ) {
			big = Number.MAX_VALUE / 4;
			return [ -big, -big ];
		}
		f = log(x) + a * log(y) + t;
		f1 = ( 1 / x ) - ( a / y );
		return [ f, f1 ];
	};
}

module.exports = temme_root_finder;
