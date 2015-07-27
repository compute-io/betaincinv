'use strict';

// MODULES //

var erfcinv = require( 'compute-erfcinv/lib/number.js' ),
	inverse_students_t_hill = require( './inverse_students_t_hill.js' ),
	inverse_students_t_body_series = require( './inverse_students_t_body_series.js' ),
	inverse_students_t_tail_series = require( './inverse_students_t_tail_series.js' ),
	ldexp = require( './ldexp.js' );


// FUNCTIONS //

var abs = Math.abs,
	acos = Math.acos,
	cos = Math.cos,
	floor = Math.floor,
	pow = Math.pow,
	round = Math.round,
	sin = Math.sin,
	sqrt = Math.sqrt;


// CONSTANTS //

var PI = Math.PI,
	ROOT_TWO = 1.414213562373095048801688724209698078e+00;


function inverse_students_t( df, u, v, pol, pexact ) {
	//
	// df = number of degrees of freedom.
	// u = probablity.
	// v = 1 - u.
	// l = lanczos type to use.
	//
	var invert, result, tmp,
		tolerance, alpha, root_alpha,
		r, x, a, b, c,
		p, p0, p2, p4, p5,
		crossover;

	invert = false;
	result = 0;
	if ( pexact ) {
		pexact.value = false;
	}
	if ( u > v ) {
		// function is symmetric, invert it:
		tmp = v;
		v = u;
		u = tmp;
		invert = true;
	}
	if( (floor(df) === df) && (df < 20) ) {
		//
		// we have integer degrees of freedom, try for the special
		// cases first:
		//
		tolerance = ldexp( 1.0, ( 2 * 53 ) / 3);

		switch( floor( df ) ) {
			case 1:
				//
				// df = 1 is the same as the Cauchy distribution, see
				// Shaw Eq 35:
				//
				if ( u === 0.5 ) {
					result = 0;
				} else {
					result = -cos( PI * u) / sin( PI * u);
				}
				if ( pexact ) {
					pexact.value = true;
				}
			break;
			case 2:
				//
				// df = 2 has an exact result, see Shaw Eq 36:
				//
				result = (2 * u - 1) / sqrt(2 * u * v);
				if ( pexact ) {
					pexact.value = true;
				}
			break;
			case 4:
				//
				// df = 4 has an exact result, see Shaw Eq 38 & 39:
				//
				alpha = 4 * u * v;
				root_alpha = sqrt( alpha );
				r = 4 * cos( acos( root_alpha ) / 3 ) / root_alpha;
				x = sqrt( r - 4 );
				result = u - 0.5 < 0 ? -x : x;
				if ( pexact ) {
					pexact.value = true;
				}
			break;
			case 6:
				//
				// We get numeric overflow in this area:
				//
				if ( u < 1e-150 ) {
					return ( invert ? -1 : 1 ) * inverse_students_t_hill(df, u, pol);
				}
				//
				// Newton-Raphson iteration of a polynomial case,
				// choice of seed value is taken from Shaw's online
				// supplement:
				//
				a = 4 * (u - u * u);//1 - 4 * (u - 0.5f) * (u - 0.5f);
				b = pow( a, 1/3 );
				c = 0.85498797333834849467655443627193;
				p = 6 * (1 + c * (1 / b - 1));
				do {
					p2 = p * p;
					p4 = p2 * p2;
					p5 = p * p4;
					p0 = p;
					// next term is given by Eq 41:
					p = 2 * (8 * a * p5 - 270 * p2 + 2187) / (5 * (4 * a * p4 - 216 * p - 243));
				} while( abs( (p - p0) / p ) > tolerance );
				//
				// Use Eq 45 to extract the result:
				//
				p = sqrt( p - df );
				result = (u - 0.5) < 0 ? -p : p;
			break;
		default:
			if ( df > 0x10000000 ) {
				result = - erfcinv( 2 * u ) * ROOT_TWO;
				if( pexact && (df >= 1e20) ) {
					pexact.value = true;
				}
			} else if ( df < 3 ) {
				//
				// Use a roughly linear scheme to choose between Shaw's
				// tail series and body series:
				//
				crossover = 0.2742 - df * 0.0242143;
				if ( u > crossover ) {
					result = inverse_students_t_body_series(df, u, pol);
				} else {
					result = inverse_students_t_tail_series(df, u, pol);
				}
			} else {
				//
				// Use Hill's method except in the exteme tails
				// where we use Shaw's tail series.
				// The crossover point is roughly exponential in -df:
				//
				crossover = ldexp( 1.0, round( df / -0.654 ) );
				if( u > crossover ) {
					result = inverse_students_t_hill(df, u, pol);
				} else {
					result = inverse_students_t_tail_series(df, u, pol);
				}
			}
		}
	} else {
		if(df > 0x10000000) {
			result = -erfcinv( 2 * u ) * ROOT_TWO;
			if ( (pexact) && (df >= 1e20) ) {
				pexact.value = true;
			}
		} else if ( df < 3 ) {
			//
			// Use a roughly linear scheme to choose between Shaw's
			// tail series and body series:
			//
			crossover = 0.2742 - df * 0.0242143;
			if ( u > crossover ) {
				result = inverse_students_t_body_series(df, u, pol);
			} else {
				result = inverse_students_t_tail_series(df, u, pol);
			}
		} else {
			//
			// Use Hill's method except in the exteme tails
			// where we use Shaw's tail series.
			// The crossover point is roughly exponential in -df:
			//
			crossover = ldexp( 1.0, round( df / -0.654) );
			if ( u > crossover ) {
				result = inverse_students_t_hill(df, u, pol);
			} else {
				result = inverse_students_t_tail_series(df, u, pol);
			}
		}
	}
	return invert ? -result : result;
}

function find_ibeta_inv_from_t_dist( a, p, py, pol ) {
	var u, v, df, t;
	u = p / 2;
	v = 1 - u;
	df = a * 2;
	t = inverse_students_t( df, u, v, pol );
	if ( py ) {
		py.value = t * t / (df + t * t);
	}
	return df / (df + t * t);
}

module.exports = find_ibeta_inv_from_t_dist;
