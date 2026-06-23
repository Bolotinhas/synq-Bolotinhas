//
// Created by adenilton on 10/4/25.
//

#include "../include/ctrl_qubit_gate.h"
#include <algorithm>

ctrl_qubit_gate::ctrl_qubit_gate(int control_qubit, int target_qubit, Eigen::Matrix2cd uMatrix) {
    control = control_qubit;
    target = target_qubit;
    data = uMatrix;
    
    decompose();
}

void ctrl_qubit_gate::decompose() {
    angles = OneQubit::zyz_decomposition(data);
}

std::string ctrl_qubit_gate::to_qasm() {
    std::string qasm;
    
    qasm += "OPENQASM 3.0;\n";
    qasm += "include \"stdgates.inc\";\n";

    int num_qubits = std::max(control, target) + 1;
    qasm += "qreg q[" + std::to_string(num_qubits) + "];\n\n";

    // Porta de fase com ângulo alfa aplicada no controle
    qasm += "p(" + std::to_string(angles.alpha) + ") q[" + std::to_string(control) + "];\n";

    // Matriz C
    qasm += "rz(" + std::to_string((angles.delta - angles.beta) / 2.0) + ") q[" + std::to_string(target) + "];\n";

    // CNOT 1
    qasm += "cx q[" + std::to_string(control) + "], q[" + std::to_string(target) + "];\n";

    // Matriz B
    qasm += "rz(" + std::to_string(-(angles.delta + angles.beta) / 2.0) + ") q[" + std::to_string(target) + "];\n";
    qasm += "ry(" + std::to_string(-angles.gamma / 2.0) + ") q[" + std::to_string(target) + "];\n";

    // CNOT 2
    qasm += "cx q[" + std::to_string(control) + "], q[" + std::to_string(target) + "];\n";

    // Matriz A
    qasm += "ry(" + std::to_string(angles.gamma / 2.0) + ") q[" + std::to_string(target) + "];\n";
    qasm += "rz(" + std::to_string(angles.beta) + ") q[" + std::to_string(target) + "];\n";

    return qasm;
}
