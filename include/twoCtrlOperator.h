/**
 * @file twoCtrlOperator.h
 * @brief Definition of the twoCtrlOperator function representing a operation controlled by 2 control bits.
 *
 * This file is part of the SynQ compiler frontend and defines the AST node
 * responsible for representing an operator controlled by two control bits applied to a state.
 *
 * @author SynQ team 
 */

#pragma once
#include <Eigen/Dense>
#include "iast.h"

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
 * @class rzNode
 * @brief AST node representing a quantum RZ rotation gate.
 *
 * This node stores a rotation angle (in radians) and an optional qubit position.
 * It is used in the AST to represent the application of an RZ gate in a quantum circuit.
 */
class CtrlOperatorNode final :public IASTnode{
public:
/**
* @brief Construct an CtrlOperator node with a given number of qubits
*
* @param n_qubits number of qubits involved in the operator
*/
explicit CtrlOperatorNode(int n_qubits);

/**
* @brief Accept a visitor according to the Visitor pattern.
* 
* @param visitor Reference to a nodeVisitor instance.
*/
void accept(nodeVisitor &visitor) override;

/**
 * @brief Retrieve the node's stored data.
 *
 * @return A return_type object containing the number of qubits.
 */
return_type get_data() override;

/**
* @brief Factory method to create a unique pointer to an CtrlOperatorNode.
*
* @param num_qubits description
* @return std::unique_ptr<IASTnode> Newly created CtrlOperator node.
*/
static std::unique_ptr<IASTnode> createCtrlOperatorNode(int num_qubits);

/**
* @brief Apply a operator controlled by one control bit on a target bit
*
* @param U Unitary matrix representing the applied operation
* @param ctrl Integer representing the position of the control bit
* @param target Integer representing the position of the target bit
*/
Eigen::MatrixXcd oneCtrlOperator(const Eigen::Matrix2cd& U, int ctrl, int target);

/**
* @brief Apply a operator controlled by two control bits on a target bit
*
* @param U Unitary matrix representing the applied operation
* @param ctrl1 Integer representing the position of the first control bit
* @param ctrl2 Integer representing the position of the second control bit
* @param target Integer representing the position of the target bit
*/
Eigen::MatrixXcd twoCtrlOperator(const Eigen::Matrix2cd& U, int ctrl1, int ctrl2, int target);

private:
/** @brief Number of qubits involved in the operation. */
int num_qubits;

};
