// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2025, the Anboto author and contributors
#include <Core/Core.h>
#include "Eigen.h"

#include <Functions4U/EnableWarnings.h>

namespace Upp {

template <>
bool MultiDimMatrix<Eigen::Matrix<double, -1, -1>>::IsNullInstance() const {
    if (IsEmpty())
		return true;
	for (int i = 0; i < size(); ++i) {
		if (d[i].size() == 0)
			return true;
		for (int j = 0; j < d[i].size(); ++j) 
			if (IsNull(d[i].array()(j)))
				return true;
	}
	return false;
}

template <>
bool MultiDimMatrix<Eigen::Matrix<double, -1, 1>>::IsNullInstance() const {
    if (IsEmpty())
		return true;
	for (int i = 0; i < size(); ++i) {
		if (d[i].size() == 0)
			return true;
		for (int j = 0; j < d[i].size(); ++j) 
			if (IsNull(d[i](j)))
				return true;
	}
	return false;
}

template <class T>
bool MultiDimMatrix<T>::IsNullInstance() const {
	if (IsEmpty())
		return true;
	for (int i = 0; i < size(); ++i) {
		if (IsNull(d[i]))
			return true;
	}
	return false;
}

}