#include<iostream>
using namespace std;
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include "../include/twoCtrlOperator.h"
#include "../include/one_qubit_gate.h"

TEST(ABCDecomposition, abc_decomposition) {
    Eigen::Matrix2cd X = OneQubit::x_matrix();
    
    // Decompõe a porta X
    abc_result abc = abc_decomposition(X);
    
    // O produto C * B * A deve ser a Identidade
    Eigen::Matrix2cd identity = Eigen::Matrix2cd::Identity();
    Eigen::Matrix2cd product = abc.C * abc.B * abc.A;

    // Print do resultado:
    std::cout << "Matrix A:\n" << abc.A << "\n\n";
    std::cout << "Matrix B:\n" << abc.B << "\n\n";
    std::cout << "Matrix C:\n" << abc.C << "\n\n";
    std::cout << "Product Matrix (C * B * A):\n" << product << "\n";

    // Verifica se o produto é numericamente igual à identidade (tolerância de 1e-6)
    EXPECT_TRUE(product.isApprox(identity, 1e-6));
}

TEST(CtrlOperators, oneCtrlOperator){
    Eigen::Matrix2cd X = OneQubit::x_matrix();
    int ctrl = 1;
    int target = 0;
    int n_qubits = 2; // Sistema de 2 qubits (dimensão 4x4)
    
    // Instanciando a classe CtrlOperatorNode passando o número de qubits
    CtrlOperatorNode operation(n_qubits);

    // Aplica CNOT 
    Eigen::MatrixXcd result = operation.oneCtrlOperator(X, ctrl, target);
    
    // A matriz esperada começa como identidade e inverte os dois últimos estados (|10> e |11>)
    Eigen::MatrixXcd expected = Eigen::MatrixXcd::Identity(4, 4);
    // Troca as linhas/colunas correspondentes aos estados 2 (|10>) e 3 (|11>)
    expected(2, 2) = 0.0; expected(2, 3) = 1.0;
    expected(3, 2) = 1.0; expected(3, 3) = 0.0;

    // Print do resultado:
    std::cout << "Result matrix:\n" << result << "\n\n";
    std::cout << "Expected matrix:\n" << expected << "\n";
    
    // Verifica se a matriz obtida é numericamente igual à matriz esperada (tolerância de 1e-6)
    EXPECT_TRUE(result.isApprox(expected, 1e-6));
}

TEST(CtrlOperators, twoCtrlOperator){
    Eigen::Matrix2cd X = OneQubit::x_matrix();
    int ctrl1 = 2;
    int ctrl2 = 1;
    int target = 0;
    int n_qubits = 3; // Sistema de 3 qubits (dimensão 8x8)
    
    // Instanciando a classe CtrlOperatorNode passando o número de qubits
    CtrlOperatorNode operation(n_qubits);

    // Executa a função de operador controlado com dois controles
    Eigen::MatrixXcd result = operation.twoCtrlOperator(X, ctrl1, ctrl2, target);
    
    // A matriz esperada começa como identidade e inverte os dois últimos estados (|110> e |111>)
    Eigen::MatrixXcd expected = Eigen::MatrixXcd::Identity(8, 8);
    // Troca as linhas/colunas correspondentes aos estados 6 (|110>) e 7 (|111>)
    expected(6, 6) = 0.0; expected(6, 7) = 1.0;
    expected(7, 6) = 1.0; expected(7, 7) = 0.0;

    // Print do resultado:
    std::cout << "Result matrix:\n" << result << "\n\n";
    std::cout << "Expected matrix:\n" << expected << "\n";
    
    // Verifica se a matriz obtida é numericamente igual à matriz esperada (tolerância de 1e-6)
    EXPECT_TRUE(result.isApprox(expected, 1e-6));
}