'use strict';

// FUNCTIONS //

var abs = Math.abs,
	pow = Math.pow,
	sign = Math.sign;

function newton_raphson_iterate( fun, guess, min, max, digits, max_iter ) {

	var result, factor, delta, delta1, delta2,
		f0 = 0,
		last_f0 = 0,
		f1,
		res, count;

	result = guess;

	factor = pow( 2, 1 - digits);
	delta = 1;
	delta1 = Number.MAX_VALUE;
	delta2 = Number.MAX_VALUE;

	count = max_iter;

	do {
		last_f0 = f0;
		delta2 = delta1;
		delta1 = delta;
		res = fun(result);
		f0 = res[ 0 ];
		f1 = res[ 1 ];
		if( 0 === f0 ) {
			break;
		}
		if( f1 === 0 ) {
			// Oops zero derivative!!!
			if( last_f0 === 0 ) {
				// this must be the first iteration, pretend that we had a
				// previous one at either min or max:
				if( result === min ) {
					guess = max;
				} else {
					guess = min;
				}
				last_f0 = fun(guess);
				delta = guess - result;
			}
			if( sign(last_f0) * sign(f0) < 0 ) {
				// we've crossed over so move in opposite direction to last step:
				if ( delta < 0 ) {
				  delta = (result - min) / 2;
				} else {
					delta = (result - max) / 2;
				}
			} else {
				// move in same direction as last step:
				if ( delta < 0 ) {
					delta = (result - max) / 2;
				} else {
					delta = (result - min) / 2;
				}
			}
	 	} else {
			delta = f0 / f1;
		}
		if( abs(delta * 2) > abs(delta2) ) {
			// last two steps haven't converged, try bisection:
			delta = (delta > 0) ? (result - min) / 2 : (result - max) / 2;
		}
		guess = result;
		result -= delta;
		if( result <= min ) {
			delta = 0.5 * (guess - min);
			result = guess - delta;
			if( (result === min) || (result === max) ) {
				break;
			}
		} else if ( result >= max ) {
			delta = 0.5 * (guess - max);
			result = guess - delta;
			if( (result === min) || (result === max) ) {
				break;
			}
		}
		// update brackets:
		if( delta > 0 ) {
			max = guess;
		} else {
			min = guess;
		}
	} while( --count && ( abs(result * factor) < abs(delta) ) );

	return result;
}

module.exports = newton_raphson_iterate;
