'use strict';

function evaluate_polynomial( poly, z, count ) {
	var sum, i;
	sum = poly[count - 1];
	for( i = count - 2; i >= 0; --i) {
		sum *= z;
		sum += poly[i];
	}
	return sum;
}

module.exports = evaluate_polynomial;
