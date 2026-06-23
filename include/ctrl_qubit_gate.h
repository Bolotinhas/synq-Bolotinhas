//
// Created by adenilton on 10/4/25.
//

/**
 * @file ctrl_qubit_gate.h
 * @brief Definition of the ctrl_qubit_gate class for testing controlled unitary decompositions.
 *
 * This file is part of the SynQ compiler frontend and defines a standalone
 * class responsible for decomposing and validating controlled unitary gates.
 *
 * @author SynQ team
 */

#ifndef TEST_TWO_QUBIT_GATE_H
#define TEST_TWO_QUBIT_GATE_H
#include <Eigen/Dense>
#include <string>
#include "one_qubit_gate.h"

/**
 * @class ctrl_qubit_gate
 * @brief Standalone class representing a controlled arbitrary unitary gate.
 *
 * This class decomposes a controlled unitary into base rotations and CNOTs 
 * generating a OpenQASM 3.0 script for simulation and testing.
 */
class ctrl_qubit_gate {
public:
    /**
     * @brief Construct a controlled unitary gate for testing.
     *
     * @param control_qubit Index of the control qubit.
     * @param target_qubit Index of the target qubit.
     * @param uMatrix The 2x2 unitary matrix to be decomposed.
     */
    explicit ctrl_qubit_gate(int control_qubit, int target_qubit, Eigen::Matrix2cd uMatrix);

    /**
     * @brief Generates a standalone OpenQASM 3.0 script representing the decomposed circuit.
     *
     * @return std::string containing the full QASM code.
     */
    std::string to_qasm();

private:
    /** @brief Index of the control qubit. */
    int control;
    
    /** @brief Index of the target qubit. */
    int target;
    
    /** @brief Original 2x2 unitary matrix. */
    Eigen::Matrix2cd data;
    
    /** @brief Extracted Euler angles from ZYZ decomposition. */
    zyz_result angles;

    /**
     * @brief Performs the ZYZ decomposition to populate the angles.
     */
    void decompose();
};


#endif //TEST_TWO_QUBIT_GATE_H