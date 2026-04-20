#ifndef AGA8_PVT_CALCULATOR_H
#define AGA8_PVT_CALCULATOR_H

#include <array>

namespace aga8_pvt {

struct GasComposition {
    std::array<double, 22> molarFractions{};  // индексы 1..21
    bool isValid(double tolerance = 1e-6) const;
};

struct ThermodynamicProperties {
    double compressibilityFactor{};
    double massDensityKgPerM3{};
    double specificEnthalpyKJPerKg{};
    double specificInternalEnergyKJPerKg;
    double specificEntropyKJPerKgK{};
    double specificCvKJPerKgK{};
    double specificCpKJPerKgK{};
    double jouleThomsonCoeffKPerMPa{};
    double adiabaticExponent{};
    double speedOfSoundMPerS{};
    double molarMassKgPerKmol{};
    double molarEnthalpyKJPerKmol{};
    void print() const;
};

class Aga8PvtCalculator {
public:
    // Расчет B_n* для заданного состава
    static std::array<double, 19> calculateBnStar(const std::array<double, 22>& composition);
    
    //Расчет второго вириального коэффициента B(τ, X)
    static double calculateB(double tau, const std::array<double, 19>& BnStar);
    
    //Расчет всех C_n для заданного состава 
    static std::array<double, 59> calculateCn(const std::array<double, 22>& composition);
    
    // Расчет смесевого параметра размера K(X)
    static double calculateK(const std::array<double, 22>& composition);
    
    // Вспомогательные методы
    static double getMolarMass(const GasComposition& composition);

    static ThermodynamicProperties calculateFromPressure(const GasComposition& composition, double pressureMpa,
        double temperatureK);
    
    static ThermodynamicProperties calculateFromDensity(
        const GasComposition& composition,
        double densityKgPerM3,
        double temperatureK);
        
    static double calculateMolarMass(const std::array<double, 22>& x);
};

} // namespace aga8_pvt

#endif