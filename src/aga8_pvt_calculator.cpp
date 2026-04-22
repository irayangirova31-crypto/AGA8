#include "aga8_pvt_calculator.h"
#include "aga8_tables.h"
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace aga8_pvt {


// ФИЗИЧЕСКИЕ КОНСТАНТЫ.

static constexpr double R_KJ = 8.314510;        // кДж/(кмоль·К).
static constexpr double R_MPa = 0.008314510;    // МПа·м³/(кмоль·К).
static constexpr double L_K = 1.0000;              // K.
static constexpr double T0_K = 298.1500;          // K.
static constexpr double P0_MPa = 0.101325;      // МПа.



// ФУНКЦИЯ 1: Расчет G_ij (уравнение D.5).

static double calculateGij(int i, int j) {
    // G_ij = G_ij* · (G_i + G_j) / 2.
    double GijStar = tables::GIJ[i][j];
    double Gi = tables::GI[i];
    double Gj = tables::GI[j];
    return GijStar * (Gi + Gj) / 2.0;
}

// ФУНКЦИЯ 2: Расчет E_ij (уравнение D.4).

static double calculateEij(int i, int j) {
    // E_ij = E_ij* · √(E_i · E_j).
    double EijStar = tables::EIJ[i][j];
    double Ei = tables::EI[i];
    double Ej = tables::EI[j];
    return EijStar * std::sqrt(Ei * Ej);
}


// ФУНКЦИЯ 3: Расчет B_nij* (уравнение D.3).

static double calculateBijStar(int n, int i, int j) {
    // 1. G_ij (уравнение D.5).
    double Gij = calculateGij(i, j);
    
    // 2. Q_i * Q_j.
    double Qij = tables::QI[i] * tables::QI[j];
    
    // 3. √(F_i · F_j).
    double Fij = std::sqrt(tables::FI[i] * tables::FI[j]);
    
    // 4. S_i · S_j.
    double Sij = tables::SI[i] * tables::SI[j];
    
    // 5. W_i · W_j.
    double Wij = tables::WI[i] * tables::WI[j];
    
    // 6. Получаем коэффициенты g_n, q_n, f_n, s_n, w_n из таблицы D.1.
    int gn = tables::GN[n];
    int qn = tables::QN[n];
    int fn = tables::FN[n];
    int sn = tables::SN[n];
    int wn = tables::WN[n];
    
    // 7. Расчет пяти членов (значение + 1 - коэффициент)^{коэффициент}.
    double term1 = (gn == 0) ? 1.0 : std::pow(Gij + 1.0 - gn, gn);
    double term2 = (qn == 0) ? 1.0 : std::pow(Qij + 1.0 - qn, qn);
    double term3 = (fn == 0) ? 1.0 : std::pow(Fij + 1.0 - fn, fn);
    double term4 = (sn == 0) ? 1.0 : std::pow(Sij + 1.0 - sn, sn);
    double term5 = (wn == 0) ? 1.0 : std::pow(Wij + 1.0 - wn, wn);
    
    return term1 * term2 * term3 * term4 * term5;
}


// ФУНКЦИЯ 3: Расчет всех B_n* (уравнение D.2).

std::array<double, 19> Aga8PvtCalculator::calculateAllBnStar(const std::array<double, 22>& x) {
    std::array<double, 19> BnStar = {0};  // индексы 1..18
    
    // Внешний цикл по n = 1..18
    for (int n = 1; n <= 18; ++n) {
        double sum = 0.0;
        double un = tables::UN[n];
        
        // Двойная сумма по компонентам i и j
        for (int i = 1; i <= 21; ++i) {
            if (x[i] == 0.0) continue;  // пропускаем компоненты с нулевой долей.
            
            for (int j = 1; j <= 21; ++j) {
                if (x[j] == 0.0) continue;  // пропускаем компоненты с нулевой долей.
                
                // 1. Расчет B_nij*
                double BijStar = calculateBijStar(n, i, j);
                // 2. Расчет E_ij (уравнение D.4) Фактор E_ij^{u_n}
                double Eij = calculateEij(i, j);
                double EijPower = std::pow(Eij, un);
                // 3. Фактор (K_i · K_j)^{3/2}
                double KiKj = tables::KI[i] * tables::KI[j];
                double Kfactor = std::pow(KiKj, 1.5);
                
                // 3. Добавляем член в сумму.
                sum += x[i] * x[j] * BijStar * EijPower * Kfactor;
            }
        }
        
        // 4. B_n* = a_n · сумма.
        BnStar[n] = tables::AN[n] * sum;
    }
    
    return BnStar;
}



// ФУНКЦИЯ: Расчет второго вириального коэффициента B(τ, X).
// Уравнение D.1: B(τ, X) = Σ_{n=1}^{18} B_n* · τ^{u_n}.

double Aga8PvtCalculator::calculateB(double tau, const std::array<double, 19>& BnStar) {
    double B = 0.0;
    
    for (int n = 1; n <= 18; ++n) {
        double un = tables::UN[n];           // степень u_n из таблицы D.1.
        double tauPow = std::pow(tau, un);   // τ^{u_n}.
        B += BnStar[n] * tauPow;
    }
    
    return B;
}


static double calculateF(const std::array<double, 22>& x) {
    double F = 0.0;
    for (int i = 1; i <= 21; ++i) {
        F += x[i] * x[i] * tables::FI[i];
    }
    return F;
}

static double calculateQ(const std::array<double, 22>& x) {
    double Q = 0.0;
    for (int i = 1; i <= 21; ++i) {
        Q += x[i] * tables::QI[i];
    }
    return Q;
}

static double calculateG(const std::array<double, 22>& x) {
    // Шаг 1: Σ x_i · G_i
    double G = 0.0;
    for (int i = 1; i <= 21; ++i) {
        G += x[i] * tables::GI[i];
    }
    
    // Шаг 2: Двойная сумма по i < j
    for (int i = 1; i <= 20; ++i) {
        if (x[i] == 0) continue;
        for (int j = i + 1; j <= 21; ++j) {
            if (x[j] == 0) continue;
            
            double xij = x[i] * x[j];
            double GijStar = tables::GIJ[i][j];
            
            G += xij * (GijStar - 1.0) * (tables::GI[i] + tables::GI[j]);
        }
    }
    
    return G;
}


static double calculateV(const std::array<double, 22>& x) {
    // Шаг 1: Σ x_i · E_i^{5/2}
    double sum1 = 0.0;
    for (int i = 1; i <= 21; ++i) {
        sum1 += x[i] * std::pow(tables::EI[i], 2.5);
    }
    
    // V⁵ = (sum1)²
    double V5 = sum1 * sum1;
    
    // Шаг 2: Двойная сумма по i < j
    for (int i = 1; i <= 20; ++i) {
        if (x[i] == 0) continue;
        for (int j = i + 1; j <= 21; ++j) {
            if (x[j] == 0) continue;
            
            double xij = x[i] * x[j];
            double Vij5 = std::pow(tables::UIJ[i][j], 5);
            double EiEj25 = std::pow(tables::EI[i] * tables::EI[j], 2.5);
            
            V5 += 2.0 * xij * (Vij5 - 1.0) * EiEj25;
        }
    }
    
    // V = (V⁵)^{1/5}
    return std::pow(V5, 0.2);
}


// РАСЧЕТ ВСЕХ C_n (уравнение D.6) для n = 13..58.

std::array<double, 59> Aga8PvtCalculator::calculateAllCn(const std::array<double, 22>& x) {
    std::array<double, 59> Cn = {0};  // индексы 1..58.
    
    // Предварительные расчеты (зависят только от состава, не от n).
    double V = calculateV(x);
    double G = calculateG(x);
    double Q = calculateQ(x);
    double F = calculateF(x);
    
    // Цикл по n = 13..58
    for (int n = 13; n <= 58; ++n) {
        int gn = tables::GN[n];
        int qn = tables::QN[n];
        int fn = tables::FN[n];
        double un = tables::UN[n];
        double an = tables::AN[n];
        

        double term1 = (gn == 0) ? 1.0 : std::pow(G + 1.0 - gn, gn);
        double term2 = (qn == 0) ? 1.0 : std::pow(Q * Q + 1.0 - qn, qn);
        double term3 = (fn == 0) ? 1.0 : std::pow(F + 1.0 - fn, fn);
        double term4 = std::pow(V, un);
        
        // C_n = a_n · term1 · term2 · term3 · term4.
        Cn[n] = an * term1 * term2 * term3 * term4;
    }
    
    return Cn;
}


// Расчет смесевого параметра размера K(X).
double Aga8PvtCalculator::calculateK(const std::array<double, 22>& x){
    // Шаг 1: Σ x_i · K_i^{5/2}.
    double sum1 = 0.0;
    for (int i = 1; i <= 21; ++i) {
        sum1 += x[i] * std::pow(tables::KI[i], 2.5);
    }
    
    // Шаг 2: K⁵ = (sum1)².
    double K5 = sum1 * sum1;
    
    // Шаг 3: Двойная сумма по i < j.
    for (int i = 1; i <= 20; ++i) {
        if (x[i] == 0) continue;
        for (int j = i + 1; j <= 21; ++j) {
            if (x[j] == 0) continue;
            
            double xij = x[i] * x[j];
            double Kij5 = std::pow(tables::KIJ[i][j], 5.0);
            double KiKj25 = std::pow(tables::KI[i] * tables::KI[j], 2.5);
            
            K5 += 2.0 * xij * (Kij5 - 1.0) * KiKj25;
        }
    }
    
    // Шаг 4: K = (K⁵)^{1/5}.
    return std::pow(K5, 0.2);
}


// ФУНКЦИЯ: Расчет фактора сжимаемости Z (уравнение 9 из ГОСТ).

static double calculateZ(
    double delta,           // относительная плотность.
    double tau,             // обратная температура (L/T).
    double B,               // второй вириальный коэффициент B(τ,X).
    double K,               // смесевой параметр размера.
    const std::array<double, 59>& Cn) {  // коэффициенты C_n.
    
    double K3 = K * K * K;
    
    // Член 1: 1 + B·δ/K³.
    double Z = 1.0 + B * delta / K3;
    
    // Член 2: -δ·Σ_{n=13}^{18} C_n·τ^{u_n}.
    double sumC13_18 = 0.0;
    for (int n = 13; n <= 18; ++n) {
        double un = tables::UN[n];
        sumC13_18 += Cn[n] * std::pow(tau, un);
    }
    Z -= delta * sumC13_18;
    
    // Член 3: Σ_{n=13}^{58} C_n·τ^{u_n}·δ^{b_n}·(b_n - c_n·k_n·δ^{k_n})·exp(-c_n·δ^{k_n}).
    for (int n = 13; n <= 58; ++n) {
        double un = tables::UN[n];
        int bn = tables::BN[n];
        int cn = tables::CN[n];
        int kn = tables::KN[n];
        
        double term = Cn[n] * std::pow(tau, un) * std::pow(delta, bn);
        double bracket = bn - cn * kn * std::pow(delta, kn);
        double expTerm = std::exp(-cn * std::pow(delta, kn));
        
        Z += term * bracket * expTerm;
    }
    
    return Z;
}


// ФУНКЦИЯ: Расчет давления по заданным δ, T, K, Z.

static double calculatePressure(
    double delta,      // относительная плотность
    double tau,        // обратная температура (L/T)
    double K,          // смесевой параметр размера
    double Z) {        // фактор сжимаемости
    
    double K3 = K * K * K;
    double L = 1.0;    // L = 1 K (константа)
    
    // p = δ · R · L / (τ · K³) · Z
    double p = delta * R_MPa * L * Z / (tau * K3);
    
    return p;
}


// РАСЧЕТ ИДЕАЛЬНО-ГАЗОВОЙ СОСТАВЛЯЮЩЕЙ φ₀ (уравнение B.3)
double Aga8PvtCalculator::calculateMolarMass(const std::array<double, 22>& x) {
    double M = 0.0;
    for (int i = 1; i <= 21; ++i) {
        M += x[i] * tables::MW[i];
    }
    return M;
}

double Aga8PvtCalculator::getMolarMass(const GasComposition& composition) {
    return calculateMolarMass(composition.molarFractions);
}
static double calculatePhi0(
    double delta,
    double tau,
    const std::array<double, 22>& x) {
    
    // Стандартные значения.
    double delta0 = P0_MPa / (R_MPa * T0_K);
    double tau0 = L_K / T0_K;
    
    // Член ln(δ/δ₀) + ln(τ₀/τ).
    double phi0 = std::log(delta / delta0) + std::log(tau0 / tau);
    
    // Сумма по компонентам.
    for (int i = 1; i <= 21; ++i) {
        if (x[i] == 0.0) continue;
        
        const auto& coeff = tables::IDEAL_COEFFS[i];
        
        double term = coeff.A01 + coeff.A02 * tau + coeff.B0 * std::log(tau);
        
        // Добавляем ln(x[i]) только если x[i] > 0.
        term += std::log(x[i]);
        
        // Члены с sinh/cosh - только если аргумент не ноль.
        if (coeff.D0 != 0.0 && coeff.D0 * tau > 0) {
            term += coeff.C0 * std::log(std::sinh(coeff.D0 * tau));
        }
        
        if (coeff.F0 != 0.0 && coeff.F0 * tau > 0) {
            term -= coeff.E0 * std::log(std::cosh(coeff.F0 * tau));
        }
        
        if (coeff.H0 != 0.0 && coeff.H0 * tau > 0) {
            term += coeff.G0 * std::log(std::sinh(coeff.H0 * tau));
        }
        
        if (coeff.J0 != 0.0 && coeff.J0 * tau > 0) {
            term -= coeff.I0 * std::log(std::cosh(coeff.J0 * tau));
        }
        
        phi0 += x[i] * term;
    }
    
    // Проверка на NaN
    if (std::isnan(phi0)) {
        return 0.0;
    }
    
    return phi0;
}

// Первая ПРОИЗВОДНАЯ φ_τ0 

static double calculatePhiTau0(double tau, const std::array<double, 22>& x) {
    double phiTau = 0.0;
    
    for (int i = 1; i <= 21; ++i) {
        if (x[i] == 0.0) continue;
        
        const auto& coeff = tables::IDEAL_COEFFS[i];
        
        double term = coeff.A02 + (coeff.B0 - 1) / tau;
        
        // Член с C0 и D0 - только если D0 != 0.
        if (coeff.D0 != 0.0) {
            term += coeff.C0 * coeff.D0 * std::cosh(coeff.D0 * tau) / std::sinh(coeff.D0 * tau);
        }
        
        // Член с E0 и F0 - только если F0 != 0;
        if (coeff.F0 != 0.0) {
            term -= coeff.E0 * coeff.F0 * std::sinh(coeff.F0 * tau) / std::cosh(coeff.F0 * tau);
        }
        
        // Член с G0 и H0 - только если H0 != 0.
        if (coeff.H0 != 0.0) {
            term += coeff.G0 * coeff.H0 * std::cosh(coeff.H0 * tau) / std::sinh(coeff.H0 * tau);
        }
        
        // Член с I0 и J0 - только если J0 != 0.
        if (coeff.J0 != 0.0) {
            term -= coeff.I0 * coeff.J0 * std::sinh(coeff.J0 * tau) / std::cosh(coeff.J0 * tau);
        }
        
        phiTau += x[i] * term;
    }
    
    return phiTau;
}


// ВТОРАЯ ПРОИЗВОДНАЯ φ_ττ0 = ∂²φ₀/∂τ² (для расчета Cv)

static double calculatePhiTauTau0(double tau, const std::array<double, 22>& x) {
    double phiTauTau = 0.0;
    
    for (int i = 1; i <= 21; ++i) {
        if (x[i] == 0.0) continue;
        
        const auto& coeff = tables::IDEAL_COEFFS[i];
        
        double term = -(coeff.B0 - 1) / (tau * tau);
        
        if (coeff.D0 != 0.0) {
            double Dtau = coeff.D0 * tau;
            term -= coeff.C0 * coeff.D0 * coeff.D0 / (std::sinh(Dtau) * std::sinh(Dtau));
        }
        
        if (coeff.F0 != 0.0) {
            double Ftau = coeff.F0 * tau;
            term -= coeff.E0 * coeff.F0 * coeff.F0 / (std::cosh(Ftau) * std::cosh(Ftau));
        }
        
        if (coeff.H0 != 0.0) {
            double Htau = coeff.H0 * tau;
            term -= coeff.G0 * coeff.H0 * coeff.H0 / (std::sinh(Htau) * std::sinh(Htau));
        }
        
        if (coeff.J0 != 0.0) {
            double Jtau = coeff.J0 * tau;
            term -= coeff.I0 * coeff.J0 * coeff.J0 / (std::cosh(Jtau) * std::cosh(Jtau));
        }
        
        phiTauTau += x[i] * term;
    }
    
    return phiTauTau;
}

// ФУНКЦИЯ: τ·φ_τ = τ·(∂φ/∂τ) (уравнение C.2)
// Используется для расчета H, S, U

static double calculateTauPhiTau(
    double delta,
    double tau,
    double K,
    const std::array<double, 59>& Cn,
    const std::array<double, 19>& BnStar,
    const std::array<double, 22>& x) {
    
    double K3 = K * K * K;
    
    // Член 1: τ·φ_0,τ (идеально-газовая часть).
    double tauPhiTau0 = tau * calculatePhiTau0(tau, x);

    // Член 2: (δ/K³)·Σ_{n=1}^{18} u_n·B_n*·τ^{u_n}.
    double sumBn = 0.0;
    for (int n = 1; n <= 18; ++n) {
        double un = tables::UN[n];
        sumBn += un * BnStar[n] * std::pow(tau, un);
    }
    
    // Член 3: -δ·Σ_{n=13}^{18} u_n·C_n·τ^{u_n}.
    double sumCn13_18 = 0.0;
    for (int n = 13; n <= 18; ++n) {
        double un = tables::UN[n];
        sumCn13_18 += un * Cn[n] * std::pow(tau, un);
    }
    
    // Член 4: Σ_{n=13}^{58} u_n·C_n·τ^{u_n}·δ^{b_n}·exp(-c_n·δ^{k_n}).
    double sumCn13_58 = 0.0;
    for (int n = 13; n <= 58; ++n) {
        double un = tables::UN[n];
        int bn = tables::BN[n];
        int cn = tables::CN[n];
        int kn = tables::KN[n];
        
        sumCn13_58 += un * Cn[n] 
                    * std::pow(tau, un)
                    * std::pow(delta, bn)
                    * std::exp(-cn * std::pow(delta, kn));
    }
    
    return tauPhiTau0 + (delta / K3) * sumBn - delta * sumCn13_18 + sumCn13_58;
}



// ФУНКЦИЯ: τ²·φ_ττ = τ²·(∂²φ/∂τ²) (уравнение C.3)
// Используется для расчета Cv


static double calculateTau2PhiTauTau(
    double delta,
    double tau,
    double K,
    const std::array<double, 59>& Cn,
    const std::array<double, 19>& BnStar,
    const std::array<double, 22>& x) {
    
    double K3 = K * K * K;
    
    // Член 1: τ²·φ_0,ττ (идеально-газовая часть).
    double tau2PhiTauTau0 = tau * tau * calculatePhiTauTau0(tau, x);
    
    // Член 2: (δ/K³)·Σ_{n=1}^{18} (u_n² - u_n)·B_n*·τ^{u_n}.
    double sumBn = 0.0;
    for (int n = 1; n <= 18; ++n) {
        double un = tables::UN[n];
        sumBn += (un * un - un) * BnStar[n] * std::pow(tau, un);
    }
    
    // Член 3: -δ·Σ_{n=13}^{18} (u_n² - u_n)·C_n·τ^{u_n}.
    double sumCn13_18 = 0.0;
    for (int n = 13; n <= 18; ++n) {
        double un = tables::UN[n];
        sumCn13_18 += (un * un - un) * Cn[n] * std::pow(tau, un);
    }
    
    // Член 4: Σ_{n=13}^{58} (u_n² - u_n)·C_n·τ^{u_n}·δ^{b_n}·exp(-c_n·δ^{k_n}).
    double sumCn13_58 = 0.0;
    for (int n = 13; n <= 58; ++n) {
        double un = tables::UN[n];
        int bn = tables::BN[n];
        int cn = tables::CN[n];
        int kn = tables::KN[n];
        
        sumCn13_58 += (un * un - un) * Cn[n] 
                    * std::pow(tau, un)
                    * std::pow(delta, bn)
                    * std::exp(-cn * std::pow(delta, kn));
    }
    
    return tau2PhiTauTau0 + (delta / K3) * sumBn - delta * sumCn13_18 + sumCn13_58;
}



// ФУНКЦИЯ: δ·φ_δ = δ·(∂φ/∂δ) (уравнение C.4).


static double calculateDeltaPhiDelta(
    double delta,
    double tau,
    double B,
    double K,
    const std::array<double, 59>& Cn) {
    
    double K3 = K * K * K;
    
    // Член 1: 1
    double result = 1.0;
    
    // Член 2: B·δ/K³
    result += B * delta / K3;
    
    // Член 3: -δ·Σ_{n=13}^{18} C_n·τ^{u_n}
    double sumCn13_18 = 0.0;
    for (int n = 13; n <= 18; ++n) {
        double un = tables::UN[n];
        sumCn13_18 += Cn[n] * std::pow(tau, un);
    }
    result -= delta * sumCn13_18;
    
    // Член 4: Σ_{n=13}^{58} C_n·τ^{u_n}·δ^{b_n}·(b_n - c_n·k_n·δ^{k_n})·exp(-c_n·δ^{k_n})
    for (int n = 13; n <= 58; ++n) {
        double un = tables::UN[n];
        int bn = tables::BN[n];
        int cn = tables::CN[n];
        int kn = tables::KN[n];
        
        double deltaPowKn = std::pow(delta, kn);
        double expTerm = std::exp(-cn * deltaPowKn);
        
        result += Cn[n] 
                * std::pow(tau, un)
                * std::pow(delta, bn)
                * (bn - cn * kn * deltaPowKn)
                * expTerm;
    }
    
    return result;
}


// ФУНКЦИЯ: Φ₁ = (δ²·φ_δ)_δ (уравнение C.5)
// Используется для расчета Cp и w

static double calculatePhi1(
    double delta,
    double tau,
    double B,
    double K,
    const std::array<double, 59>& Cn) {
    
    double K3 = K * K * K;
    
    // Член 1: 1
    double result = 1.0;
    
    // Член 2: 2·B·δ/K³
    result += 2.0 * B * delta / K3;
    
    // Член 3: -2·δ·Σ_{n=13}^{18} C_n·τ^{u_n}
    double sumCn13_18 = 0.0;
    for (int n = 13; n <= 18; ++n) {
        double un = tables::UN[n];
        sumCn13_18 += Cn[n] * std::pow(tau, un);
    }
    result -= 2.0 * delta * sumCn13_18;
    
    // Член 4: Σ C_n·τ^{u_n}·δ^{b_n}·[b_n - (1+k_n)·c_n·k_n·δ^{k_n} + (b_n - c_n·k_n·δ^{k_n})²]·exp(...)
    for (int n = 13; n <= 58; ++n) {
        double un = tables::UN[n];
        int bn = tables::BN[n];
        int cn = tables::CN[n];
        int kn = tables::KN[n];
        
        double deltaPowKn = std::pow(delta, kn);
        double expTerm = std::exp(-cn * deltaPowKn);
        
        double bracket = bn - (1.0 + kn) * cn * kn * deltaPowKn
                       + std::pow(bn - cn * kn * deltaPowKn, 2);
        
        result += Cn[n] 
                * std::pow(tau, un)
                * std::pow(delta, bn)
                * bracket
                * expTerm;
    }
    
    return result;
}


// ФУНКЦИЯ: Φ₂ = -τ²·(δ·φ_δ/τ)_τ (уравнение C.6).
// Используется для расчета Cp и μ.


static double calculatePhi2(
    double delta,
    double tau,
    double K,
    const std::array<double, 59>& Cn,
    const std::array<double, 19>& BnStar) {
    
    double K3 = K * K * K;
    
    // Член 1: 1.
    double result = 1.0;
    
    // Член 2: (δ/K³)·Σ_{n=1}^{18} (1-u_n)·B_n*·τ^{u_n}.
    double sumBn = 0.0;
    for (int n = 1; n <= 18; ++n) {
        double un = tables::UN[n];
        sumBn += (1.0 - un) * BnStar[n] * std::pow(tau, un);
    }
    result += (delta / K3) * sumBn;
    
    // Член 3: -δ·Σ_{n=13}^{18} (1-u_n)·C_n·τ^{u_n}.
    double sumCn13_18 = 0.0;
    for (int n = 13; n <= 18; ++n) {
        double un = tables::UN[n];
        sumCn13_18 += (1.0 - un) * Cn[n] * std::pow(tau, un);
    }
    result -= delta * sumCn13_18;
    
    // Член 4: Σ_{n=13}^{58} (1-u_n)·C_n·τ^{u_n}·δ^{b_n}·(b_n - c_n·k_n·δ^{k_n})·exp(-c_n·δ^{k_n}).
    for (int n = 13; n <= 58; ++n) {
        double un = tables::UN[n];
        int bn = tables::BN[n];
        int cn = tables::CN[n];
        int kn = tables::KN[n];
        
        double deltaPowKn = std::pow(delta, kn);
        double expTerm = std::exp(-cn * deltaPowKn);
        
        result += (1.0 - un) * Cn[n] 
                * std::pow(tau, un)
                * std::pow(delta, bn)
                * (bn - cn * kn * deltaPowKn)
                * expTerm;
    }
    
    return result;
}


// РАСЧЕТ НЕИДЕАЛЬНОЙ СОСТАВЛЯЮЩЕЙ φᵣ (уравнение 11)

static double calculatePhiR(
    double delta,
    double tau,
    double B,
    double K,
    const std::array<double, 59>& Cn) {
    
    double K3 = K * K * K;
    
    // Член 1: B·δ/K³.
    double phiR = B * delta / K3;
    
    // Член 2: -δ·Σ_{n=13}^{18} C_n·τ^{u_n}.
    double sumC13_18 = 0.0;
    for (int n = 13; n <= 18; ++n) {
        double un = tables::UN[n];
        sumC13_18 += Cn[n] * std::pow(tau, un);
    }
    phiR -= delta * sumC13_18;
    
    // Член 3: Σ_{n=13}^{58} C_n·τ^{u_n}·δ^{b_n}·exp(-c_n·δ^{k_n}).
    for (int n = 13; n <= 58; ++n) {
        double un = tables::UN[n];
        int bn = tables::BN[n];
        int cn = tables::CN[n];
        int kn = tables::KN[n];
        
        double term = Cn[n] 
                    * std::pow(tau, un) 
                    * std::pow(delta, bn)
                    * std::exp(-cn * std::pow(delta, kn));
        
        phiR += term;
    }
    
    return phiR;
}




// ЕДИНАЯ ФУНКЦИЯ: расчет всех термодинамических свойств по δ и T.
// Используется как для calculateFromPressure, так и для calculateFromDensity.

static ThermodynamicProperties calculateAllProperties(
    const std::array<double, 22>& x,
    double delta,
    double temperatureK,
    double K,
    const std::array<double, 59>& Cn,
    const std::array<double, 19>& BnStar,
    double pressureMpa = -1.0) {  // -1 = давление неизвестно.
    
    double tau = 1.0 / temperatureK;
    double M = Aga8PvtCalculator::calculateMolarMass(x);
    double K3 = K * K * K;
    
    // 1. Базовые величины.
    double B = Aga8PvtCalculator::calculateB(tau, BnStar);
    
    // 2. Расчет производных по формулам C.2 - C.6.
    double tauPhiTau = calculateTauPhiTau(delta, tau, K, Cn, BnStar, x);
    double tau2PhiTauTau = calculateTau2PhiTauTau(delta, tau, K, Cn, BnStar, x);
    double deltaPhiDelta = calculateDeltaPhiDelta(delta, tau, B, K, Cn);
    double Phi1 = calculatePhi1(delta, tau, B, K, Cn);
    double Phi2 = calculatePhi2(delta, tau, K, Cn, BnStar);
    
    // 3. Фактор сжимаемости Z (уравнение 17).
    double Z = deltaPhiDelta;
    
    // 4. Если давление не задано, вычисляем его.
    double p = pressureMpa;
    if (p < 0) {
        p = delta * R_MPa * L_K * Z / (tau * K3);
    }
    
    // 5. Плотность.
    double rhoMol = delta / K3;
    double rhoMass = rhoMol * M;
    
    // 6. Полная свободная энергия φ (для энтропии).
    double phi0 = calculatePhi0(delta, tau, x);
    double phiR = calculatePhiR(delta, tau, B, K, Cn);
    double phi = phi0 + phiR;


    // 7. Термодинамические свойства.
    double uMolar = R_KJ * temperatureK * tauPhiTau;                             // (19)
    double hMolar = R_KJ * temperatureK * (tauPhiTau + deltaPhiDelta);           // (20)
    double sMolar = R_KJ * (tauPhiTau - phi);                                    // (21)
    double cvMolar = -R_KJ * tau2PhiTauTau;                                      // (22)
    double cpMolar = cvMolar + R_KJ * Phi2 * Phi2 / Phi1;                        // (23)
    double kappa = (Phi1 / Z) * (cpMolar / cvMolar);                             // (25)
    double mu = (Phi2 / Phi1 - 1.0) / (cpMolar * (p / (Z * R_KJ * temperatureK)));// (24)
    double w = std::sqrt(Z * kappa * R_KJ * 1000.0 * temperatureK / M);          // (26)
    
    // 8. Заполнение структуры.
    ThermodynamicProperties props;
    props.compressibilityFactor = Z;
    props.massDensity = rhoMass;
    props.specificInternalEnergy = uMolar / M;
    props.specificEnthalpy = hMolar / M;
    props.specificEntropy = sMolar / M;
    props.specificCv = cvMolar / M;
    props.specificCp = cpMolar / M;
    props.jouleThomsonCoeff = mu;
    props.adiabaticExponent = kappa;
    props.speedOfSound = w;
    props.molarMass = M;
    props.molarEnthalpy = hMolar;
    
    return props;
}

void ThermodynamicProperties::print() const {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Термодинамические свойства газа" << std::endl;
    std::cout << "  Z = " << compressibilityFactor << std::endl;
    std::cout << "  D = " << massDensity << " кг/м³" << std::endl;
    std::cout << "  U = " << specificInternalEnergy << " кДж/кг" << std::endl;
    std::cout << "  H = " << specificEnthalpy << " кДж/кг" << std::endl;
    std::cout << "  S = " << specificEntropy << " кДж/(кг·К)" << std::endl;
    std::cout << "  cv = " << specificCv << " кДж/(кг·К)" << std::endl;
    std::cout << "  cp = " << specificCp << " кДж/(кг·К)" << std::endl;
    std::cout << "  μ = " << jouleThomsonCoeff << " К/МПа" << std::endl;
    std::cout << "  κ = " << adiabaticExponent << std::endl;
    std::cout << "  w = " << speedOfSound << " м/с" << std::endl;
    std::cout << "  M = " << molarMass << " кг/кмоль" << std::endl;
}

// ФУНКЦИЯ: Итерационное решение для δ по заданному давлению.

static double solveForDeltaFromPressure(
    const std::array<double, 22>& x,
    double pressureMpa, // целевое давление (известно).
    double temperatureK,
    double K,
    const std::array<double, 59>& Cn,
    const std::array<double, 19>& BnStar) {
    
    double tau = 1.0 / temperatureK;
    double B = Aga8PvtCalculator::calculateB(tau, BnStar);
    
    // Начальное приближение δ (идеальный газ: Z=1).
    double K3 = K * K * K;
    double delta = pressureMpa * tau * K3 / (R_MPa * 1.0);
    
    const double tolerance = 1e-8;   // точность по давлению, МПа.
    const int maxIter = 50;
    
    
    for (int iter = 0; iter < maxIter; ++iter) {
        // 1. Рассчитываем Z по текущему δ.
        double Z = calculateZ(delta, tau, B, K, Cn);
        
        // 2. Рассчитываем давление по текущему δ и Z.
        double pCalc = calculatePressure(delta, tau, K, Z);
        
        // 3. Вычисляем ошибку.
        double error = pCalc - pressureMpa;
        
        // 4. Проверка сходимости.
        if (std::abs(error) < tolerance) {
            break;
        }
        
        // 5. Коррекция δ (метод простой итерации).
        delta = delta * pressureMpa / pCalc;
    }
    return delta;
}



bool GasComposition::isValid(double tolerance) const {
    double sum = 0.0;
    for (int i = 1; i <= 21; ++i) {
        sum += molarFractions[i];
    }
    return std::abs(sum - 1.0) < tolerance;
}

// Расчет по давлению и температуре

ThermodynamicProperties Aga8PvtCalculator::calculateFromPressure(
    const GasComposition& composition,
    double pressureMpa,
    double temperatureK) {
    
    // 1. Проверка входных данных.
    if (!composition.isValid()) {
        throw std::invalid_argument("Sum of molar fractions must be 1.0");
    }
    if (pressureMpa <= 0.0 || pressureMpa > 30.0) {
        throw std::out_of_range("Pressure must be in (0, 30] MPa");
    }
    if (temperatureK < 250.0 || temperatureK > 350.0) {
        throw std::out_of_range("Temperature must be in [250, 350] K");
    }
    
    // 2. Рассчитываем коэффициенты.
    auto BnStar = Aga8PvtCalculator::calculateAllBnStar(composition.molarFractions);
    auto Cn = Aga8PvtCalculator::calculateAllCn(composition.molarFractions);
    double K = calculateK(composition.molarFractions);
    
    // 3. Итерационное решение для δ
    double delta = solveForDeltaFromPressure(
        composition.molarFractions,
        pressureMpa,
        temperatureK,
        K, Cn, BnStar);
    
    // Единый расчет свойств
    return calculateAllProperties(
        composition.molarFractions, delta, temperatureK,
        K, Cn, BnStar, pressureMpa);
}

ThermodynamicProperties Aga8PvtCalculator::calculateFromDensity(
    const GasComposition &composition,
     double densityKgPerM3,
     double temperatureK){
    
    // 1. Проверка входных данных.
    if (!composition.isValid()) {
        throw std::invalid_argument("Sum of molar fractions must be 1.0");
    }
    if (densityKgPerM3 <= 0.0) {
        throw std::out_of_range("Density must be positive");
    }
    if (temperatureK < 250.0 || temperatureK > 350.0) {
        throw std::out_of_range("Temperature must be in [250, 350] K");
    }
    

    auto BnStar = Aga8PvtCalculator::calculateAllBnStar(composition.molarFractions);
    auto Cn = Aga8PvtCalculator::calculateAllCn(composition.molarFractions);
    double K = Aga8PvtCalculator::calculateK(composition.molarFractions);
    double M = getMolarMass(composition);
    
    // 3. ПРЯМОЙ РАСЧЕТ δ.
    //δ = D · K³ / M
    double K3 = K * K * K;
    double delta = densityKgPerM3 * K3 / M;
    
    // Единый расчет свойств (давление будет вычислено внутри).
    return calculateAllProperties(
        composition.molarFractions, delta, temperatureK,
        K, Cn, BnStar, -1.0);  // -1 = давление неизвестно.
}

} // namespace aga8_pvt