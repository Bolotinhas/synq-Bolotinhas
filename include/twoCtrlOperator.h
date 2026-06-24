/**
 * @file twoCtrlOperator.h
 * @brief Definition of the twoCtrlOperator function representing a operation controlled by 2 control bits.
 *
 * This file is part of the SynQ compiler frontend and defines the AST node
 * responsible for representing an RZ rotation applied to a specific qubit.
 *
 * @author SynQ team 
 */

#pragma once
#include <Eigen/Dense>

/**
 * @struct abc_result
 * @brief Struct used to store ABC matrix decomposition and delta angle 
 */
struct abc_result {
    Eigen::Matrix2cd A;
    Eigen::Matrix2cd B;
    Eigen::Matrix2cd C;
    double delta;
};

/**
* @brief Construct an ABC matrix decomposition on a unitary matrix V based on Barenco Theorem
*
* @param V Unitary matrix that is decomposed 
*/
abc_result abc_decomposition(const Eigen::Matrix2cd& V);

/**
* @brief Apply a operator controlled by one control bit on a target bit
*
* @param U Unitary matrix representing the applied operation
* @param ctrl Integer representing the position of the control bit
* @param target Integer representing the position of the target bit
* @param n_qubits Number of qubits involved in the operation
*/
Eigen::MatrixXcd oneCtrlOperator(const Eigen::Matrix2cd& U, int ctrl, int target, int n_qubits);

/**
* @brief Apply a operator controlled by two control bits on a target bit
*
* @param U Unitary matrix representing the applied operation
* @param ctrl1 Integer representing the position of the first control bit
* @param ctrl2 Integer representing the position of the second control bit
* @param target Integer representing the position of the target bit
* @param n_qubits Number of qubits involved in the operation
*/
Eigen::MatrixXcd twoCtrlOperator(const Eigen::Matrix2cd& U, int ctrl1, int ctrl2, int target, int n_qubits);
