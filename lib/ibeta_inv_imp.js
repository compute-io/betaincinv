'use strict';

// MODULES //

var beta = require( 'compute-beta' ),
	betainc = require( 'compute-betainc/lib/ibeta_imp.js' ),
	polyval = require( 'compute-polynomial' ),
	halley_iterate = require( './halley_iterate.js' ),
	ibeta_roots = require( './ibeta_roots.js' ),
	temme_method_1_ibeta_inverse = require( './temme1.js' ),
	temme_method_2_ibeta_inverse = require( './temme2.js' ),
	temme_method_3_ibeta_inverse = require( './temme3.js' ),
	find_ibeta_inv_from_t_dist = require( './find_ibeta_inv_from_t_dist.js' );


// FUNCTIONS //

var log1p = Math.log1p || function(x) {
  return Math.log(1 + x);
};

var expm1 = Math.expm1 || function(x) {
	return Math.exp(x) - 1;
};

var abs = Math.abs,
	asin = Math.asin,
	exp = Math.exp,
	log = Math.log,
	min = Math.min,
	max = Math.max,
	pow = Math.pow,
	sin = Math.sin,
	sqrt = Math.sqrt;


// CONSTANTS //

var HALF_PI = Math.PI / 2,
	EPSILON = Number.EPSILON || 2.2204460492503130808472633361816e-16;


//  Copyright John Maddock 2006.
//  Copyright Paul A. Bristow 2007
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

function ibeta_inv_imp( a, b, p, q ) {
	var tmp, lx, x, y, yp,
		invert,
		lower, upper, l, u,
		minv, maxv,
		r, theta, lambda,
		ps, fs,
		ppa, bet, xs, xs2, xg,
		digits, max_iter;
	//
	// The flag invert is set to true if we swap a for b and p for q,
	// in which case the result has to be subtracted from 1:
	//
	invert = false;
	//
	// Handle trivial cases first:
	//
	if( q === 0 ) {
		return [ 1, 0 ];
	} else if( p === 0 ) {
		return [ 0, 1 ];
	} else if( a === 1 ) {
		if( b === 1 ) {
			return [ p, 1 - p ];
		}
		// Change things around so we can handle as b == 1 special case below:
		tmp = b;
		b = a;
		a = tmp;

		tmp = q;
		q = p;
		p = tmp;

		invert = true;
	}
	//
	// Depending upon which approximation method we use, we may end up
	// calculating either x or y initially (where y = 1-x):
	//
	x = 0; // Set to a safe zero to avoid a
	// MSVC 2005 warning C4701: potentially uninitialized local variable 'x' used
	// But code inspection appears to ensure that x IS assigned whatever the code path.

	// For some of the methods we can put tighter bounds
	// on the result than simply [0,1]:
	//
	lower = 0;
	upper = 1;
	//
	// Student's T with b = 0.5 gets handled as a special case, swap
	// around if the arguments are in the "wrong" order:
	//
	if ( a === 0.5 ) {
		if ( b === 0.5 ) {
			x = sin( p * HALF_PI );
			x *= x;
			y = sin( q * HALF_PI );
			y *= y;
			return [ x, y ];
		} else if( b > 0.5 ) {

			tmp = b;
			b = a;
			a = tmp;

			tmp = q;
			q = p;
			p = tmp;

			invert = !invert;
		}
	}
	//
	// Select calculation method for the initial estimate:
	//
	if( (b === 0.5) && (a >= 0.5) && (p !== 1) ) {
		//
		// We have a Student's T distribution:
		yp = {};
		x = find_ibeta_inv_from_t_dist( a, p, yp );
		y = yp.value;
	} else if ( b === 1 ) {
		if ( p < q ) {
			if( a > 1 ) {
				x = pow(p, 1 / a);
				y = -expm1(log(p) / a );
			} else {
				x = pow(p, 1 / a);
				y = 1 - x;
			}
		} else {
			x = exp( log1p(-q) / a);
			y = -expm1( log1p(-q) / a );
		}
		if ( invert ) {
			tmp = y;
			y = x;
			x = tmp;
		}
		return [ x, y ];
	} else if ( a + b > 5 ) {
		//
		// When a+b is large then we can use one of Prof Temme's
		// asymptotic expansions, begin by swapping things around
		// so that p < 0.5, we do this to avoid cancellations errors
		// when p is large.
		//
		if( p > 0.5 ) {
			tmp = b;
			b = a;
			a = tmp;

			tmp = q;
			q = p;
			p = tmp;

			invert = !invert;
		}
		minv = min(a, b);
		maxv = max(a, b);
		if( (sqrt(minv) > (maxv - minv)) && (minv > 5) ) {
			//
			// When a and b differ by a small amount
			// the curve is quite symmetrical and we can use an error
			// function to approximate the inverse. This is the cheapest
			// of the three Temme expantions, and the calculated value
			// for x will never be much larger than p, so we don't have
			// to worry about cancellation as long as p is small.
			//
			x = temme_method_1_ibeta_inverse( a, b, p );
			y = 1 - x;
		} else {
			r = a + b;
			theta = asin(sqrt(a / r));
			lambda = minv / r;
			if( (lambda >= 0.2) && (lambda <= 0.8) && (r >= 10) ) {
				//
				// The second error function case is the next cheapest
				// to use, it brakes down when the result is likely to be
				// very small, if a+b is also small, but we can use a
				// cheaper expansion there in any case.  As before x won't
				// be much larger than p, so as long as p is small we should
				// be free of cancellation error.
				//
				ppa = pow(p, 1/a);
				if( (ppa < 0.0025) && (a + b < 200) ) {
					x = ppa * pow(a * beta(a, b ), 1/a);
				} else {
					x = temme_method_2_ibeta_inverse( p, r, theta );
				}
				y = 1 - x;
		} else {
			//
			// If we get here then a and b are very different in magnitude
			// and we need to use the third of Temme's methods which
			// involves inverting the incomplete gamma.  This is much more
			// expensive than the other methods.  We also can only use this
			// method when a > b, which can lead to cancellation errors
			// if we really want y (as we will when x is close to 1), so
			// a different expansion is used in that case.
			//
			if ( a < b ) {
				tmp = b;
				b = a;
				a = tmp;

				tmp = q;
				q = p;
				p = tmp;
				invert = !invert;
			}
			//
			// Try and compute the easy way first:
			//
			bet = 0;
			if ( b < 2 ) {
				bet = beta( a, b );
			}
			if ( bet !== 0 ) {
				y = pow(b * q * bet, 1/b);
				x = 1 - y;
			} else {
				y = 1;
			}
			if( y > 1e-5 ) {
				x = temme_method_3_ibeta_inverse( a, b, p, q );
				y = 1 - x;
			}
		}
	 }
	} else if((a < 1) && (b < 1)) {
		//
		// Both a and b less than 1,
		// there is a point of inflection at xs:
		//
		xs = (1 - a) / (2 - a - b);
		//
		// Now we need to ensure that we start our iteration from the
		// right side of the inflection point:
		//
		fs = betainc( a, b, xs ) - p;
		if ( abs(fs) / p < EPSILON * 3 ) {
			// The result is at the point of inflection, best just return it:
			if ( invert ) {
				return [ 1 - xs, xs ];
			} else {
				return [ xs, 1 - xs ];
			}
		}
		if ( fs < 0 ) {
			tmp = b;
			b = a;
			a = tmp;

			tmp = q;
			q = p;
			p = tmp;

			invert = !invert;
			xs = 1 - xs;
		}
		xg = pow( a * p * beta( a, b ), 1/a );
		x = xg / ( 1 + xg );
		y = 1 / ( 1 + xg );
		//
		// And finally we know that our result is below the inflection
		// point, so set an upper limit on our search:
		//
		if ( x > xs ) {
			x = xs;
		}
		upper = xs;
	} else if( (a > 1) && (b > 1) ) {
		//
		// Small a and b, both greater than 1,
		// there is a point of inflection at xs,
		// and it's complement is xs2, we must always
		// start our iteration from the right side of the
		// point of inflection.
		//
		xs = (a - 1) / (a + b - 2);
		xs2 = (b - 1) / (a + b - 2);
		ps = betainc( a, b, xs ) - p;

		if ( ps < 0 ) {
			tmp = b;
			b = a;
			a = tmp;

			tmp = q;
			q = p;
			p = tmp;

			tmp = xs2;
			xs2 = xs;
			xs = tmp;

			invert = !invert;
		}
		//
		// Estimate x and y, using expm1 to get a good estimate
		// for y when it's very small:
		//
		lx = log(p * a * beta( a, b ) ) / a;
		x = exp( lx );
		y = x < 0.9 ? 1 - x : -expm1(lx);

		if( (b < a) && (x < 0.2) ) {
			//
			// Under a limited range of circumstances we can improve
			// our estimate for x, frankly it's clear if this has much effect!
			//
			var ap1, bm1, a_2, a_3, b_2, terms;

			ap1 = a - 1;
			bm1 = b - 1;
			a_2 = a * a;
			a_3 = a * a_2;
			b_2 = b * b;
			terms = new Array( 5 );
			terms[0] = 0;
			terms[1] = 1;
			terms[2] = bm1 / ap1;
			ap1 *= ap1;
			terms[3] = bm1 * (3 * a * b + 5 * b + a_2 - a - 4) / (2 * (a + 2) * ap1);
			ap1 *= (a + 1);
			terms[4] = bm1 * (33 * a * b_2 + 31 * b_2 + 8 * a_2 * b_2 - 30 * a * b - 47 * b + 11 * a_2 * b + 6 * a_3 * b + 18 + 4 * a - a_3 + a_2 * a_2 - 10 * a_2);
			terms[4] /= (3 * (a + 3) * (a + 2) * ap1);
			x = polyval( terms, x );
		}
		//
		// And finally we know that our result is below the inflection
		// point, so set an upper limit on our search:
		//
		if ( x > xs ) {
			x = xs;
		}
		upper = xs;
	} else /*if((a <= 1) != (b <= 1))*/ {
		//
		// If all else fails we get here, only one of a and b
		// is above 1, and a+b is small.  Start by swapping
		// things around so that we have a concave curve with b > a
		// and no points of inflection in [0,1].  As long as we expect
		// x to be small then we can use the simple (and cheap) power
		// term to estimate x, but when we expect x to be large then
		// this greatly underestimates x and leaves us trying to
		// iterate "round the corner" which may take almost forever...
		//
		// We could use Temme's inverse gamma function case in that case,
		// this works really rather well (albeit expensively) even though
		// strictly speaking we're outside it's defined range.
		//
		// However it's expensive to compute, and an alternative approach
		// which models the curve as a distorted quarter circle is much
		// cheaper to compute, and still keeps the number of iterations
		// required down to a reasonable level.  With thanks to Prof Temme
		// for this suggestion.
		//
		if(b < a) {
			tmp = b;
			b = a;
			a = tmp;

			tmp = q;
			q = p;
			p = tmp;
			invert = !invert;
		}
		if( pow(p, 1/a) < 0.5 ) {
			x = pow(p * a * beta( a, b ), 1 / a);
			if ( x === 0 ) {
				x = Number.MIN_VALUE;
			}
		 	y = 1 - x;
		}
		else /*if(pow(q, 1/b) < 0.1)*/ {
			// model a distorted quarter circle:
			y = pow( 1 - pow( p, b * beta( a, b ) ), 1/b );
			if(y === 0) {
				y = Number.MIN_VALUE;
			}
			x = 1 - y;
		}
	}
	//
	// Now we have a guess for x (and for y) we can set things up for
	// iteration.  If x > 0.5 it pays to swap things round:
	//
	if ( x > 0.5 ) {
		tmp = b;
		b = a;
		a = tmp;

		tmp = q;
		q = p;
		p = tmp;

		tmp = y;
		y = x;
		x = tmp;

		invert = !invert;
		l = 1 - upper;
		u = 1 - lower;
		lower = l;
		upper = u;
	}
	//
	// lower bound for our search:
	//
	// We're not interested in denormalised answers as these tend to
	// these tend to take up lots of iterations, given that we can't get
	// accurate derivatives in this area (they tend to be infinite).
	//
	if( lower === 0 ) {
		if( invert ) {
			//
			// We're not interested in answers smaller than machine epsilon:
			//
			lower = EPSILON;
			if ( x < lower ) {
				x = lower;
			}
		} else {
			lower = Number.MIN_VALUE;
		}
		if ( x < lower ) {
			x = lower;
		}
	}
	//
	// Figure out how many digits to iterate towards:
	//
	digits = 100 / 2;
	if( (x < 1e-50) && ( (a < 1) || (b < 1) ) ) {
		//
		// If we're in a region where the first derivative is very
		// large, then we have to take care that the root-finder
		// doesn't terminate prematurely.  We'll bump the precision
		// up to avoid this, but we have to take care not to set the
		// precision too high or the last few iterations will just
		// thrash around and convergence may be slow in this case.
		// Try 3/4 of machine epsilon:
		//
		digits *= 3;
		digits /= 2;
	}
	//
	// Now iterate, we can use either p or q as the target here
	// depending on which is smaller:
	//
	max_iter = 1000;
	x = halley_iterate(
		ibeta_roots(a, b, (p < q ? p : q), (p < q ? false : true)), x, lower, upper, digits, max_iter );
	// Tidy up, if we "lower" was too high then zero is the best answer we have:
	//
	if( x === lower ) {
		x = 0;
	}
	if ( invert ) {
		return [ 1-x, x ];
	}
	return [ x, 1-x ];
}

module.exports= ibeta_inv_imp;
