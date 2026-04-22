#ifndef AGA8_PVT_CALCULATOR_H
#define AGA8_PVT_CALCULATOR_H

#include <array>

namespace aga8_pvt {

// СОСТАВ ГАЗА

// Индексы массива соответствуют следующему порядку компонентов:
//
// Индекс | Компонент          | Хим. формула
// -------|--------------------|-------------
//   1    | Азот               | N2
//   2    | Диоксид углерода   | CO2
//   3    | Метан              | CH4
//   4    | Этан               | C2H6
//   5    | Пропан             | C3H8
//   6    | н-Бутан            | n-C4H10
//   7    | Изобутан           | i-C4H10
//   8    | н-Пентан           | n-C5H12
//   9    | Изопентан          | i-C5H12
//  10    | н-Гексан           | n-C6H14
//  11    | н-Гептан           | n-C7H16
//  12    | н-Октан            | n-C8H18
//  13    | н-Нонан            | n-C9H20
//  14    | н-Декан            | n-C10H22
//  15    | Водород            | H2
//  16    | Кислород           | O2
//  17    | Оксид углерода     | CO
//  18    | Вода               | H2O
//  19    | Сероводород        | H2S
//  20    | Гелий              | He
//  21    | Аргон              | Ar
//
// Пример использования:
//   GasComposition comp;
//   comp.molarFractions[3] = 0.95;   // 95% метана
//   comp.molarFractions[4] = 0.03;   // 3% этана
//   comp.molarFractions[1] = 0.02;   // 2% азота
//   // Сумма всех долей должна быть равна 1.0
//

struct GasComposition {
    std::array<double, 22> molarFractions{};  // индексы 1..21.
    bool isValid(double tolerance = 1e-6) const;
};

struct ThermodynamicProperties {
    double compressibilityFactor{};     // Z - фактор сжимаемости
    double massDensity{};               // D - массовая плотность, кг/м³
    double specificEnthalpy{};          // H - удельная энтальпия, кДж/кг
    double specificInternalEnergy;      // U - удельная внутренняя энергия, кДж/кг
    double specificEntropy{};           // S - удельная энтропия, кДж/(кг·К)
    double specificCv{};                // cv - удельная изохорная теплоемкость, кДж/(кг·К)
    double specificCp{};                // cp - удельная изобарная теплоемкость, кДж/(кг·К)
    double jouleThomsonCoeff{};         // μ - коэффициент Джоуля-Томсона, К/МПа
    double adiabaticExponent{};         // κ - показатель адиабаты
    double speedOfSound{};              // w - скорость звука, м/с
    double molarMass{};                 // M - молярная масса, кг/кмоль
    double molarEnthalpy{};             // h - молярная энтальпия, кДж/кмоль
    void print() const;
};

class Aga8PvtCalculator {
public:
    // Расчет B_n* для заданного состава.
    static std::array<double, 19> calculateAllBnStar(const std::array<double, 22>& composition);
    
    //Расчет второго вириального коэффициента B(τ, X).
    static double calculateB(double tau, const std::array<double, 19>& BnStar);
    
    //Расчет всех C_n для заданного состава. 
    static std::array<double, 59> calculateAllCn(const std::array<double, 22>& composition);
    
    // Расчет смесевого параметра размера K(X).
    static double calculateK(const std::array<double, 22>& composition);
    
    // Вспомогательные методы.
    static double getMolarMass(const GasComposition& composition);
    static double calculateMolarMass(const std::array<double, 22>& x);

    static ThermodynamicProperties calculateFromPressure(const GasComposition& composition, double pressureMpa,
        double temperatureK);
    
    static ThermodynamicProperties calculateFromDensity(
        const GasComposition& composition,
        double densityKgPerM3,
        double temperatureK);
        
};

} // namespace aga8_pvt.

#endif