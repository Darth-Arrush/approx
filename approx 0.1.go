/*package approx contains methods to approximate roots of polynomials*/
/*Use Taylor/Maclaurin series to approximate roots of non-polynomials*/
/*v 0.1*/

package approx

import (
	"errors"
	"math"
)

/*Implementation of polynomial, along with execute(poly Polynomial, x float64)*/
/*Enter as c₀, c₁, c₂, c₃, c₄ where cₙ is coefficient of xⁿ*/

type Polynomial []float64

func execute(poly Polynomial, x float64) float64 {
	result := 0.0
	for i := 0; i < len(poly); i += 1 {
		result += poly[i] * (math.Pow(x, float64(i)))
	}
	return result
}

/*Methods*/

/*Bracketing*/
func Bisection(poly Polynomial, x0 float64, x1 float64, accuracy float64) (float64, error) /*Finds root to y = 0 ± accuracy*/ {
	if math.Abs(execute(poly, (x0+x1)/2)) < math.Abs(accuracy) {
		return (x0 + x1) / 2, nil
	}
	if execute(poly, x0) < 0 && execute(poly, x1) > 0 {
		for math.Abs(execute(poly, (x0+x1)/2)) >= math.Abs(accuracy) {
			if execute(poly, (x0+x1)/2) > 0 {
				x1 = (x0 + x1) / 2
			} else if execute(poly, (x0+x1)/2) < 0 {
				x0 = (x0 + x1) / 2
			} else {
				return 0.0, errors.New("error")
			}
		}
		return (x0 + x1) / 2, nil
	} else if execute(poly, x0) > 0 && execute(poly, x1) < 0 {
		for math.Abs(execute(poly, (x0+x1)/2)) >= math.Abs(accuracy) {
			if execute(poly, (x0+x1)/2) > 0 {
				x0 = (x0 + x1) / 2
			} else if execute(poly, (x0+x1)/2) < 0 {
				x1 = (x0 + x1) / 2
			} else {
				return 0.0, errors.New("error")
			}
		}
	}
	return 0.0, errors.New("both x0 and x1 have same sign")
}

func BisectionIter(poly Polynomial, x0 float64, x1 float64, iterations int) (float64, error) /*Finds root in iterations iterations*/ {
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

func RegulaFalsi(poly Polynomial, x0 float64, x1 float64, accuracy float64) (float64, error) /*Finds root to y = 0 ± accuracy*/

func RegulaFalsiIter(poly Polynomial, x0 float64, x1 float64, iterations int) (float64, error) /*Finds root in iterations iterations*/

/*Interpolation*/
func Secant(poly Polynomial, x0 float64, x1 float64, accuracy float64) float64 /*Finds root to y = 0 ± accuracy*/ {
	approximations := []float64{x0, x1}
	for math.Abs(approximations[len(approximations)-1]) < math.Abs(accuracy) {
		newApproximation := ((approximations[len(approximations)-2] * execute(poly, approximations[len(approximations)-1])) - (approximations[len(approximations)-1] * execute(poly, approximations[len(approximations)-2]))) / (execute(poly, approximations[len(approximations)-1]) - execute(poly, approximations[len(approximations)-2]))
		approximations = append(approximations, newApproximation)
	}
	return approximations[len(approximations)-1]
}

func SecantIter(poly Polynomial, x0 float64, x1 float64, iterations int) float64 /*Finds root in iterations iterations*/ {
	approximations := []float64{x0, x1}
	for i := 2; i < iterations+2; i += 1 {
		newApproximation := ((approximations[i-2] * execute(poly, approximations[i-1])) - (approximations[i-1] * execute(poly, approximations[i-2]))) / (execute(poly, approximations[i-1]) - execute(poly, approximations[i-2]))
		approximations = append(approximations, newApproximation)
	}
	return approximations[len(approximations)-1]
}

/*Iterative*/
func Newton(poly Polynomial, x0 float64, x1 float64, accuracy float64) float64 /*Finds root to y = 0 ± accuracy*/ {
	var derivative Polynomial = []float64{}
	for i := 1; i < len(poly); i += 1 {
		derivative = append(derivative, float64(i)*poly[i])
	}
	approximations := []float64{poly[0]}
	for math.Abs(execute(poly, approximations[len(approximations)-1])) < math.Abs(accuracy) { /*Uses formula x_(n+1) = x_n - f(x_n)/f'(x_n)*/
		newApproximation := approximations[len(approximations)-1] - execute(poly, approximations[len(approximations)-1])/execute(derivative, approximations[len(approximations)-1])
		approximations = append(approximations, newApproximation)
	}
	return approximations[len(approximations)-1]
}

func NewtonIter(poly Polynomial, x0 float64, x1 float64, iterations int) float64 /*Finds root in iterations iterations*/ {
	var derivative Polynomial = []float64{}
	for i := 1; i < len(poly); i += 1 {
		derivative = append(derivative, float64(i)*poly[i])
	}
	approximations := []float64{poly[len(poly)-1]}
	for j := 0; j <= iterations; j += 1 { /*Uses formula x_(n+1) = x_n - f(x_n)/f'(x_n)*/
		newApproximation := approximations[j] - execute(poly, approximations[len(approximations)-1])/execute(derivative, approximations[len(approximations)-1])
		approximations = append(approximations, newApproximation)
	}
	return approximations[len(approximations)-1]
}

func Steffensen(poly Polynomial, x0 float64, x1 float64, h float64, accuracy float64) float64 /*Finds root to y = 0 ± accuracy*/ /*Redundant ∵ same as Newton, except with h -/-> 0*/

func SteffensenIter(poly Polynomial, x0 float64, x1 float64, h float64, iterations int) float64 /*Finds root in iterations iterations*/ /*Redundant ∵ same as Newton, except with h -/-> 0*/

func Halley(poly Polynomial, x0 float64, x1 float64, accuracy float64) float64 /*Finds root to y = 0 ± accuracy*/

func HalleyIter(poly Polynomial, x0 float64, x1 float64, iterations int) float64 /*Finds root in iterations iterations*/

func Broyden(poly Polynomial, x0 float64, x1 float64, accuracy float64) float64 /*Finds root to y = 0 ± accuracy*/

func BroydenIter(poly Polynomial, x0 float64, x1 float64, iterations int) float64 /*Finds root in iterations iterations*/
