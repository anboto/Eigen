// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2025, the Anboto author and contributors
#include <Core/Core.h>
#include "Eigen.h"


namespace Upp {

using namespace Eigen;

bool NonLinearOptimization(VectorXd &y, Eigen::Index numData, 
			Function <int(const VectorXd &b, VectorXd &residual)> Residual,
			double xtol, double ftol, int maxfev, int &ret) {
	BasicNonLinearOptimizationFunctor functor(Residual);
	functor.unknowns = y.size();
	functor.datasetLen = numData;
	Eigen::NumericalDiff<BasicNonLinearOptimizationFunctor> numDiff(functor);
	Eigen::LevenbergMarquardt<Eigen::NumericalDiff<BasicNonLinearOptimizationFunctor>> lm(numDiff);
	if (!IsNull(xtol))
		lm.parameters.xtol *= xtol;
	if (!IsNull(ftol))
		lm.parameters.ftol *= ftol;
	if (!IsNull(maxfev))
		lm.parameters.maxfev = maxfev;
	ret = lm.minimize(y);
	if (ret == Eigen::LevenbergMarquardtSpace::ImproperInputParameters || 
		ret == Eigen::LevenbergMarquardtSpace::TooManyFunctionEvaluation ||
		ret == Eigen::LevenbergMarquardtSpace::CosinusTooSmall) 
		return false;
	return true;
}

bool NonLinearOptimization(VectorXd &y, Eigen::Index numData, 
			Function <int(const VectorXd &b, VectorXd &residual)> Residual,
			double xtol, double ftol, int maxfev) {
	int ret;
	return NonLinearOptimization(y, numData, Residual, xtol, ftol, maxfev, ret);
}

String NonLinearOptimizationError(int error) {
	switch(error) {
	case Eigen::LevenbergMarquardtSpace::NotStarted:						return "StatusNotStarted";
	case Eigen::LevenbergMarquardtSpace::Running:							return "Running";
	case Eigen::LevenbergMarquardtSpace::ImproperInputParameters:			return "Improper input parameters";
	case Eigen::LevenbergMarquardtSpace::RelativeReductionTooSmall:			return "Relative reduction too small";
	case Eigen::LevenbergMarquardtSpace::RelativeErrorTooSmall:				return "Relative error too small";
	case Eigen::LevenbergMarquardtSpace::RelativeErrorAndReductionTooSmall:	return "Relative error and reduction too small";
	case Eigen::LevenbergMarquardtSpace::CosinusTooSmall:					return "Cosinus too small";
	case Eigen::LevenbergMarquardtSpace::TooManyFunctionEvaluation:			return "Too many function evaluations";
	case Eigen::LevenbergMarquardtSpace::FtolTooSmall:						return "Ftol too small";
	case Eigen::LevenbergMarquardtSpace::XtolTooSmall:						return "Xtol too small";
	case Eigen::LevenbergMarquardtSpace::GtolTooSmall:						return "Gtol too small";
	case Eigen::LevenbergMarquardtSpace::UserAsked:							return "User asked";
	}
	return "Unknown code";
}

bool SolveNonLinearEquations(VectorXd &y, Function <int(const VectorXd &b, VectorXd &residual)> Residual,
			double xtol, int maxfev, double factor) {
	BasicNonLinearOptimizationFunctor functor(Residual);
	HybridNonLinearSolver<BasicNonLinearOptimizationFunctor> solver(functor);
	if (!IsNull(xtol))
		solver.parameters.xtol *= xtol;
	if (!IsNull(maxfev))
		solver.parameters.maxfev = maxfev;
	if (!IsNull(factor))
		solver.parameters.factor = factor;
	int ret = solver.solveNumericalDiff(y);
	if (ret == HybridNonLinearSolverSpace::ImproperInputParameters ||
	    ret == HybridNonLinearSolverSpace::TooManyFunctionEvaluation ||
	    ret == HybridNonLinearSolverSpace::NotMakingProgressJacobian ||
	    ret == HybridNonLinearSolverSpace::NotMakingProgressIterations)
		return false;
	return true;
}

double SolveNonLinearEquation(double y, Function <double(double b)> Residual, double xtol, int maxfev, double factor) {
	VectorXd x(1), res;
	x[0] = y;
	
	auto Residual2 = [&](const VectorXd &b, VectorXd &residual)->int {
		residual[0] = Residual(b[0]);
       	if (IsNull(residual[0]))
       		return 1;
       	return 0;
	};
	if (SolveNonLinearEquations(x, Residual2, xtol, maxfev, factor))
		return x[0];

	return Null;		
}

Matrix3d SkewSymmetricMatrix(const Vector3d& r) {
    Matrix3d s;
    s << 0, -r.z(), r.y(),
         r.z(), 0, -r.x(),
        -r.y(), r.x(), 0;
    return s;
}

}