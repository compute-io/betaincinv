'use strict';

// MODULES //

var polyval = require( 'compute-polynomial' ),
	tgamma_delta_ratio = require( './tgamma_delta_ratio.js' );


// FUNCTIONS //

var pow = Math.pow,
	sqrt = Math.sqrt;

// CONSTANTS //

var PI = Math.PI;

//
// Tail and body series are due to Shaw:
//
// www.mth.kcl.ac.uk/~shaww/web_page/papers/Tdistribution06.pdf
//
// Shaw, W.T., 2006, "Sampling Student's T distribution - use of
// the inverse cumulative distribution function."
// Journal of Computational Finance, Vol 9 Issue 4, pp 37-73, Summer 2006
//
function inverse_students_t_tail_series( df, v, pol ) {
	var w, np2, np4, np6,
		rn, div, power, result;
	// Tail series expansion, see section 6 of Shaw's paper.
	// w is calculated using Eq 60:
	w = tgamma_delta_ratio(df / 2, 0.5, pol) * sqrt( df * PI ) * v;
	// define some variables:
	np2 = df + 2;
	np4 = df + 4;
	np6 = df + 6;
	//
	// Calculate the coefficients d(k), these depend only on the
	// number of degrees of freedom df, so at least in theory
	// we could tabulate these for fixed df, see p15 of Shaw:
	//
	var d = new Array( 7 );
	d[ 0 ] = 1;
	d[ 1 ] = -(df + 1) / (2 * np2);
	np2 *= (df + 2);
	d[ 2 ] = -df * (df + 1) * (df + 3) / (8 * np2 * np4);
	np2 *= df + 2;
	d[ 3 ] = -df * (df + 1) * (df + 5) * (((3 * df) + 7) * df -2) / (48 * np2 * np4 * np6);
	np2 *= (df + 2);
	np4 *= (df + 4);
	d[ 4 ] = -df * (df + 1) * (df + 7) *
		( (((((15 * df) + 154) * df + 465) * df + 286) * df - 336) * df + 64 ) /
		(384 * np2 * np4 * np6 * (df + 8));
	np2 *= (df + 2);
	d[ 5 ] = -df * (df + 1) * (df + 3) * (df + 9) *
		(((((((35 * df + 452) * df + 1573) * df + 600) * df - 2020) * df) + 928) * df -128) /
		(1280 * np2 * np4 * np6 * (df + 8) * (df + 10));
	np2 *= (df + 2);
	np4 *= (df + 4);
	np6 *= (df + 6);
	d[ 6 ] = -df * (df + 1) * (df + 11) *
		((((((((((((945 * df) + 31506) * df + 425858) * df + 2980236) * df + 11266745) * df + 20675018) * df + 7747124) * df - 22574632) * df - 8565600) * df + 18108416) * df - 7099392) * df + 884736) /
		(46080 * np2 * np4 * np6 * (df + 8) * (df + 10) * (df +12));
	//
	// Now bring everthing together to provide the result,
	//
	// this is Eq 62 of Shaw:
	rn = sqrt( df );
	div = pow( rn * w, 1 / df );
	power = div * div;
	result = polyval( d, power );
	result *= rn;
	result /= div;
	return -result;
}

module.exports = inverse_students_t_tail_series;
