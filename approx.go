/* package approx contains methods to approximate roots of real-valued polynomials */
/* Use Taylor/Maclaurin series to approximate roots of non-polynomials */
/* approx uses float64, so slight inaccuracy is to be expected */
/* v 1.0.0 */

package approx

import (
	"errors"
	"math"
)

/* Type */
type Polynomial []float64 /* Implementation of polynomial */ /* Enter as c₀, c₁, c₂, c₃, c₄, ... */

/* Auxiliary Functions */
func execute(poly Polynomial, x float64) float64 { /* Returns float64 output of poly Polynomial at x float64 */
	result := 0.0
	for i := 0; i < len(poly); i += 1 {
		result += poly[i] * (math.Pow(x, float64(i)))
	}
	return result
}

func g(poly Polynomial, x float64) float64 { /* For Steffensen and SteffensenIter */
	h := execute(poly, x)
	var result float64 = (execute(poly, (x+h)) - execute(poly, x)) / h
	return result
}

/* Methods */
/* Bracketing */
func Bisection(poly Polynomial, x0 float64, x1 float64, accuracy float64) (float64, error) /* Finds root to y = 0 ± accuracy */ {
	if math.Abs(execute(poly, (x0+x1)/2)) < math.Abs(accuracy) {
		return (x0 + x1) / 2, nil
	}
	if execute(poly, x0) < 0 && execute(poly, x1) > 0 {
		if execute(poly, (x0+x1)/2) > 0 {
			x1 = (x0 + x1) / 2
		} else if execute(poly, (x0+x1)/2) < 0 {
			x0 = (x0 + x1) / 2
		} else {
			return 0.0, errors.New("error")
		}
		Bisection(poly, x0, x1, accuracy)
		return (x0 + x1) / 2, nil
	} else if execute(poly, x0) > 0 && execute(poly, x1) < 0 {
		if execute(poly, (x0+x1)/2) > 0 {
			x0 = (x0 + x1) / 2
		} else if execute(poly, (x0+x1)/2) < 0 {
			x1 = (x0 + x1) / 2
		} else {
			return 0.0, errors.New("error")
		}
		Bisection(poly, x0, x1, accuracy)
		return (x0 + x1) / 2, nil
	}
	return 0.0, errors.New("both x0 and x1 have same sign")
}

func BisectionIter(poly Polynomial, x0 float64, x1 float64, iterations int) (float64, error) /* Finds root in iterations iterations */ {
	if execute(poly, x0) < 0 && execute(poly, x1) > 0 {
		for i := 0; i < iterations; i += 1 {
			if execute(poly, (x0+x1)/2) > 0 {
				x1 = (x0 + x1) / 2
			} else if execute(poly, (x0+x1)/2) < 0 {
				x0 = (x0 + x1) / 2
			} else if execute(poly, (x0+x1)/2) == 0 {
				break
			} else {
				return 0.0, errors.New("error")
			}
		}
		return (x0 + x1) / 2, nil
	} else if execute(poly, x0) > 0 && execute(poly, x1) < 0 {
		for i := 0; i < iterations; i += 1 {
			if execute(poly, (x0+x1)/2) > 0 {
				x0 = (x0 + x1) / 2
			} else if execute(poly, (x0+x1)/2) < 0 {
				x1 = (x0 + x1) / 2
			} else if execute(poly, (x0+x1)/2) == 0 {
				break
			} else {
				return 0.0, errors.New("error")
			}
		}
	}
	return 0.0, errors.New("both x0 and x1 have same sign")
}

/* Iterative */
func Newton(poly Polynomial, x0 float64, accuracy float64) (float64, error) /* Finds root to y = 0 ± accuracy */ {
	if len(poly) <= 1 {
		return 0.0, errors.New("enter a differentiable polynomial")
	}
	var derivative Polynomial = []float64{}
	for i := 1; i < len(poly); i += 1 {
		derivative = append(derivative, float64(i)*poly[i])
	}
	approximations := []float64{x0}
	for math.Abs(execute(poly, approximations[len(approximations)-1])) > math.Abs(accuracy) { /* Uses formula x_(n+1) = x_n - f(x_n)/f'(x_n) */
		newApproximation := approximations[len(approximations)-1] - execute(poly, approximations[len(approximations)-1])/execute(derivative, approximations[len(approximations)-1])
		approximations = append(approximations, newApproximation)
	}
	return approximations[len(approximations)-1], nil
}

func NewtonIter(poly Polynomial, x0 float64, iterations int) (float64, error) /* Finds root in iterations iterations */ {
	if len(poly) <= 1 {
		return 0.0, errors.New("enter a differentiable polynomial")
	}
	var derivative Polynomial = []float64{}
	for i := 1; i < len(poly); i += 1 {
		derivative = append(derivative, float64(i)*poly[i])
	}
	approximations := []float64{x0}
	for j := 0; j <= iterations; j += 1 { /* Uses formula x_(n+1) = x_n - f(x_n)/f'(x_n) */
		newApproximation := approximations[j] - execute(poly, approximations[len(approximations)-1])/execute(derivative, approximations[len(approximations)-1])
		approximations = append(approximations, newApproximation)
	}
	return approximations[len(approximations)-1], nil
}

func Halley(poly Polynomial, x0 float64, accuracy float64) (float64, error) /* Finds root to y = 0 ± accuracy */ {
	if len(poly) <= 2 {
		return 0.0, errors.New("enter a twice-differentiable polynomial")
	}
	var derivative Polynomial = []float64{}
	for i := 1; i < len(poly); i += 1 {
		derivative = append(derivative, float64(i)*poly[i])
	}
	var doubleDerivative Polynomial = []float64{}
	for j := 1; j < len(derivative); j += 1 {
		doubleDerivative = append(doubleDerivative, float64(j)*derivative[j])
	}
	approximations := []float64{x0}
	for math.Abs(execute(poly, approximations[len(approximations)-1])) > math.Abs(accuracy) {
		newApproximation := approximations[len(approximations)-1] - 2*execute(poly, approximations[len(approximations)-1])*execute(derivative, approximations[len(approximations)-1])/((2*execute(derivative, approximations[len(approximations)-1])*execute(derivative, approximations[len(approximations)-1]))-execute(poly, approximations[len(approximations)-1])*execute(doubleDerivative, approximations[len(approximations)-1]))
		approximations = append(approximations, newApproximation)
	}
	return approximations[len(approximations)-1], nil
}

func HalleyIter(poly Polynomial, x0 float64, iterations int) (float64, error) /* Finds root in iterations iterations */ {
	if len(poly) <= 2 {
		return 0.0, errors.New("enter a twice-differentiable polynomial")
	}
	var derivative Polynomial = []float64{}
	for i := 1; i < len(poly); i += 1 {
		derivative = append(derivative, float64(i)*poly[i])
	}
	var doubleDerivative Polynomial = []float64{}
	for j := 1; j < len(derivative); j += 1 {
		doubleDerivative = append(doubleDerivative, float64(j)*derivative[j])
	}
	approximations := []float64{x0}
	for k := 0; k < iterations; k += 1 {
		newApproximation := approximations[len(approximations)-1] - 2*execute(poly, approximations[len(approximations)-1])*execute(derivative, approximations[len(approximations)-1])/((2*execute(derivative, approximations[len(approximations)-1])*execute(derivative, approximations[len(approximations)-1]))-execute(poly, approximations[len(approximations)-1])*execute(doubleDerivative, approximations[len(approximations)-1]))
		approximations = append(approximations, newApproximation)
	}
	return approximations[len(approximations)-1], nil
}

func Steffensen(poly Polynomial, x0 float64, accuracy float64) float64 /* Finds root to y = 0 ± accuracy */ {
	approximations := []float64{x0}
	for math.Abs(execute(poly, approximations[len(approximations)-1])) > math.Abs(accuracy) { /* Uses formula x_(n+1) = x_n - f(x_n)/g(x_n) */
		newApproximation := approximations[len(approximations)-1] - execute(poly, approximations[len(approximations)-1])/g(poly, approximations[len(approximations)-1])
		approximations = append(approximations, newApproximation)
	}
	return approximations[len(approximations)-1]
}

func SteffensenIter(poly Polynomial, x0 float64, iterations int) float64 /* Finds root in iterations iterations */ {
	/*VERY slow convergence if x0 not close to root*/
	approximations := []float64{x0}
	for i := 0; i < iterations; i += 1 {
		newApproximation := approximations[len(approximations)-1] - execute(poly, approximations[len(approximations)-1])/g(poly, approximations[len(approximations)-1])
		approximations = append(approximations, newApproximation)
	}
	return approximations[len(approximations)-1]
}
