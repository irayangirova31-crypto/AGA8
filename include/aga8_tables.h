#ifndef AGA8_TABLES_H
#define AGA8_TABLES_H

#include <array>

namespace aga8_pvt {
namespace tables {


// Таблица D.1 - Коэффициенты уравнения состояния (n = 1..58).
extern const std::array<double, 59> AN;   // a_n
extern const std::array<int, 59> BN;      // b_n
extern const std::array<int, 59> CN;      // c_n
extern const std::array<int, 59> KN;      // k_n
extern const std::array<double, 59> UN;   // u_n
extern const std::array<int, 59> GN;      // g_n
extern const std::array<int, 59> QN;      // q_n
extern const std::array<int, 59> FN;      // f_n
extern const std::array<int, 59> SN;      // s_n
extern const std::array<int, 59> WN;      // w_n



// Таблица D.2 - Параметры компонентов (21 компонент).

extern const std::array<double, 22> EI;   // энергетический параметр E_i
extern const std::array<double, 22> KI;   // параметр размера K_i
extern const std::array<double, 22> GI;   // ориентационный параметр G_i
extern const std::array<double, 22> QI;   // квадрупольный параметр Q_i
extern const std::array<double, 22> FI;   // высокотемпературный параметр F_i
extern const std::array<double, 22> SI;   // дипольный параметр S_i
extern const std::array<double, 22> WI;   // параметр ассоциации W_i
extern const std::array<double, 22> MW;   // молярные массы


// Таблица B.1 - Идеально-газовые коэффициенты (21 компонент).

struct IdealGasCoeffs {
    double A01;
    double A02;
    double B0;
    double C0;
    double D0;
    double E0;
    double F0;
    double G0;
    double H0;
    double I0;
    double J0;
};

extern const std::array<IdealGasCoeffs, 22> IDEAL_COEFFS;


// Таблица D.3 - Параметры бинарного взаимодействия (21x21).

extern const std::array<std::array<double, 22>, 22> EIJ;   // E_ij*
extern const std::array<std::array<double, 22>, 22> UIJ;   // U_ij (V_ij)
extern const std::array<std::array<double, 22>, 22> KIJ;   // K_ij
extern const std::array<std::array<double, 22>, 22> GIJ;   // G_ij*

} // namespace tables
} // namespace aga8_pvt

#endif
