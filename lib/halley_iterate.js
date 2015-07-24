'use strict';

// FUNCTIONS //

var abs = Math.abs,
	max = Math.max,
	pow = Math.pow,
	sign = Math.sign;

// CONSTANTS //

var MAX_VALUE = Number.MAX_VALUE;

function halley_iterate( fun, guess, minimum, maximum, digits, max_iter ) {
	var f0 = 0, f1, f2,
		last_f0,
		result, factor,
		delta, delta1, delta2,
		out_of_bounds_sentry = false,
		res,
		count,
		num, denom,
		convergence, diff;

	result = guess;

	factor =  pow( 2, 1 - digits);
	delta = max( 10000000 * guess, 10000000 );  // arbitarily large delta
	last_f0 = 0;
	delta1 = delta;
	delta2 = delta;

	count = max_iter;

	do {
		last_f0 = f0;
		delta2 = delta1;
		delta1 = delta;
		res = fun( result);
		f0 = res[ 0 ];
		f1 = res[ 1 ];
		f2 = res[ 2 ];

		if ( f0 === 0 ) {
			break;
		}
		if (f1 === 0) {
			// Oops zero derivative!!!
			if( last_f0 === 0 ) {
				// this must be the first iteration, pretend that we had a
				// previous one at either min or max:
				if( result === minimum ) {
					guess = maximum;
				} else {
					guess = minimum;
				}
				last_f0 = fun(guess);
				delta = guess - result;
			}
			if( sign(last_f0) * sign(f0) < 0 ) {
				// we've crossed over so move in opposite direction to last step:
				if ( delta < 0 ) {
				  delta = (result - minimum) / 2;
				} else {
					delta = (result - maximum) / 2;
				}
			} else {
				// move in same direction as last step:
				if ( delta < 0 ) {
					delta = (result - maximum) / 2;
				} else {
					delta = (result - minimum) / 2;
				}
			}
		} else {
			if( f2 !== 0 ) {
				denom = 2 * f0;
				num = 2 * f1 - f0 * (f2 / f1);
				if ( (abs(num) < 1) && ( abs(denom) >= abs(num) * MAX_VALUE ) ) {
					// possible overflow, use Newton step:
					delta = f0 / f1;
				} else {
					delta = denom / num;
				}
				if( delta * f1 / f0 < 0 ) {
					// Oh dear, we have a problem as Newton and Halley steps
					// disagree about which way we should move.  Probably
					// there is cancelation error in the calculation of the
					// Halley step, or else the derivatives are so small
					// that their values are basically trash.  We will move
					// in the direction indicated by a Newton step, but
					// by no more than twice the current guess value, otherwise
					// we can jump way out of bounds if we're not careful.
					// See https://svn.boost.org/trac/boost/ticket/8314.
					delta = f0 / f1;
					if( abs(delta) > 2 * abs(guess) ) {
						delta = (delta < 0 ? -1 : 1) * 2 * abs(guess);
					}
				}
			} else {
				delta = f0 / f1;
			}
		}
		convergence = abs(delta / delta2);
		if ( (convergence > 0.8) && (convergence < 2) ) {
			// last two steps haven't converged, try bisection:
			delta = (delta > 0) ? (result - minimum) / 2 : (result - maximum) / 2;
			if ( abs(delta) > result ) {
				delta = sign(delta) * result; // protect against huge jumps!
			}
			// reset delta2 so that this branch will *not* be taken on the
			// next iteration:
			delta2 = delta * 3;
		}
		guess = result;
		result -= delta;

		// check for out of bounds step:
		if (result < minimum) {
			diff = ( (abs(minimum) < 1) && (abs(result) > 1) && ( MAX_VALUE / abs(result) < abs(minimum)) ) ? 1000  : result / minimum;
			if ( abs(diff) < 1 ) {
				diff = 1 / diff;
			}
			if ( !out_of_bounds_sentry && (diff > 0) && (diff < 3) ) {
				// Only a small out of bounds step, lets assume that the result
				// is probably approximately at minimum:
				delta = 0.99 * (guess  - minimum);
				result = guess - delta;
				out_of_bounds_sentry = true; // only take this branch once!
			} else {
				delta = (guess - minimum) / 2;
				result = guess - delta;
				if( (result === minimum) || (result === maximum) ) {
					break;
				}
			}
		} else if (result > maximum) {
			diff = ( (abs(maximum) < 1) && (abs(result) > 1) && ( MAX_VALUE / abs(result) < abs(maximum) ) ) ? 1000 : result / maximum;
			if ( abs(diff) < 1 ) {
				diff = 1 / diff;
			}
			if ( !out_of_bounds_sentry && (diff > 0) && (diff < 3) ) {
				// Only a small out of bounds step, lets assume that the result
				// is probably approximately at minimum:
				delta = 0.99 * (guess  - maximum);
				result = guess - delta;
				out_of_bounds_sentry = true; // only take this branch once!
			} else {
				delta = (guess - maximum) / 2;
				result = guess - delta;
				if( (result === minimum) || (result === maximum) ) {
					break;
				}
			}
		}
		// update brackets:
		if (delta > 0) {
			maximum = guess;
		} else {
			minimum = guess;
		}
	} while ( --count && ( abs(result * factor) < abs(delta) ) );

	return result;
}

module.exports = halley_iterate;
