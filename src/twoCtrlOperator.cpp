#include "../include/twoCtrlOperator.h"
#include "../include/one_qubit_gate.h"
#include <cmath>
#include <unsupported/Eigen/MatrixFunctions>
#include "../include/nodeVisitor.h"
using namespace std;

// Faz a decomposição de um operador unitário nas matrizes ABC
abc_result abc_decomposition(const Eigen::Matrix2cd& V){
    // Faz a decomposição zyz
    zyz_result zyz = OneQubit::zyz_decomposition(V);
    double alpha = zyz.alpha;
    double beta  = zyz.beta;
    double gamma = zyz.gamma;
    double delta = zyz.delta;

    // Constrói as matrizes A, B, C 
    Eigen::Matrix2cd A = OneQubit::rz_matrix(alpha) * OneQubit::ry_matrix(beta / 2.0);
    Eigen::Matrix2cd B = OneQubit::ry_matrix(-beta / 2.0) * OneQubit::rz_matrix(-(alpha + gamma) / 2.0);
    Eigen::Matrix2cd C = OneQubit::rz_matrix((gamma - alpha) / 2.0);

    // Monta resultado
    abc_result result;
    result.A = A;
    result.B = B;
    result.C = C;
    result.delta = delta;

    return result;
}

// Inicializa a variável privada 'num_qubits' com o valor recebido por parâmetro
CtrlOperatorNode::CtrlOperatorNode(int n_qubits) : num_qubits(n_qubits) {}

// Faz o visitor marcar esse node como visitado
void CtrlOperatorNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

// Retorna o número de qubits de um CtrlOperatorNode
return_type CtrlOperatorNode::get_data() {
    return static_cast<double>(num_qubits);
}

// Cria um ponteiro único para esse node
std::unique_ptr<IASTnode> CtrlOperatorNode::createCtrlOperatorNode(int num_qubits) {
    return std::make_unique<CtrlOperatorNode>(num_qubits);
}

// Aplica uma porta U de 1 qubit no target, controlada por um qubit ctrl
Eigen::MatrixXcd CtrlOperatorNode::oneCtrlOperator(const Eigen::Matrix2cd& U, int ctrl, int target) {
    int dim = pow(2, num_qubits); // 2^num_qubits
    Eigen::MatrixXcd result = Eigen::MatrixXcd::Identity(dim, dim); //retorna identidade para caso o controle nao for setado

    // Varre todas as colunas 
    for (int col = 0; col < dim; col++)
    {
        // Extrai apenas o bit do qubit de controle
        // Exemplo: col = 1101 e ctrl = 2 (da esquerda para direita)
        // 1101 -> 0011, após isso faz o and bit a bit
        bool c = (col >> ctrl) & 1;

        // Se o bit de controle não for 1, pula 
        if (!c) continue;

        // Identifica os índices dos estados irmãos (onde o target é |0> e |1>)
        // Serve para alterar apenas o bit target enquanto mantém os outros bits
        // base0 terá exatamente os mesmos bits de col exceto no bit de target em que é forçado a ser 0
        // base1 terá exatamente os mesmos bits de col exceto no bit de target em que é forçado a ser 1
        int base0 = col & ~(1 << target);
        int base1 = base0 | (1 << target);

        // Mapeia a matriz U de 2x2 para dentro da matriz global
        for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            int row = (i == 0) ? base0 : base1;
            int c_idx = (j == 0) ? base0 : base1;

            result(row, c_idx) = U(i, j);
        }
    }

    return result;
}

/*
Operador controlado por 2 bits de controle onde V² = U:

  q1 (ctrl1)    ───●───       ───────────●───────────●───●───
                   │                     │           │   │
  q2 (ctrl2)    ───●───   =   ───●───────⊕──●───────⊕──│───
                   │             │           │           │
  q3 (target)   ──[U]──       ──[V]─────────[V†]────────[V]──
*/

// Aplica uma porta U de 1 qubit no target, controlada por dois qubit ctrl1 e ctrl2
Eigen::MatrixXcd CtrlOperatorNode::twoCtrlOperator(const Eigen::Matrix2cd& U, int ctrl1, int ctrl2, int target) {
    // V é a raiz quadrada de U (V * V = U)
    Eigen::Matrix2cd V = U.sqrt();
    // V† é a conjugada transposta de V
    Eigen::Matrix2cd V_dagger = V.adjoint();
    // Operador NOT
    Eigen::Matrix2cd X = OneQubit::x_matrix();

    // V controlada por ctrl2 no target
    Eigen::MatrixXcd op1 = oneCtrlOperator(V, ctrl2, target);
    
    // CNOT controlada por ctrl1 no ctrl2
    Eigen::MatrixXcd op2 = oneCtrlOperator(X, ctrl1, ctrl2);
    
    // V† controlada por ctrl2 no target
    Eigen::MatrixXcd op3 = oneCtrlOperator(V_dagger, ctrl2, target);
    
    // CNOT controlada por ctrl1 no ctrl2
    Eigen::MatrixXcd op4 = oneCtrlOperator(X, ctrl1, ctrl2);
    
    // V controlada por ctrl1 no target
    Eigen::MatrixXcd op5 = oneCtrlOperator(V, ctrl1, target);

    return op5 * op4 * op3 * op2 * op1;
}
