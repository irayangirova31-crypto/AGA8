#include "aga8_pvt_calculator.h"
#include "aga8_tables.h"
#include <iostream>
#include <iomanip>

using namespace aga8_pvt;

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "  AGA8 PVT Calculator" << std::endl;
    std::cout << "  ГОСТ Р 8.662-2009" << std::endl;
    std::cout << "========================================\n" << std::endl;
 
// СОСТАВ ГАЗА

// Индексы массива соответствуют следующему порядку компонентов:
//
// Индекс | Компонент          | Хим. формула
// --------|-------------------|-------------
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



    // Состав 1 природного газа.
    GasComposition comp1;
    comp1.molarFractions[1] = 0.003000;   // Азот
    comp1.molarFractions[2] = 0.006;   // Диоксид углерода
    comp1.molarFractions[3] = 0.965;   // Метан
    comp1.molarFractions[4] = 0.018;   // Этан
    comp1.molarFractions[5] = 0.0045;   // Пропан
    comp1.molarFractions[6] = 0.001;   // н-Бутан
    comp1.molarFractions[7] = 0.001;   // Изопутан
    comp1.molarFractions[8] = 0.0003;   // н-Пентан
    comp1.molarFractions[9] = 0.0005;   // Изопентан
    comp1.molarFractions[10] = 0.0007;   // н-Гексан


    // Состав 6 природного газа.
    GasComposition comp6;
    comp6.molarFractions[1] = 0.117266;   // Азот
    comp6.molarFractions[2] = 0.011093;   // Диоксид углерода
    comp6.molarFractions[3] = 0.825198;   // Метан
    comp6.molarFractions[4] = 0.034611;   // Этан
    comp6.molarFractions[5] = 0.007645;   // Пропан
    comp6.molarFractions[6] = 0.002539;   // н-Бутан
    comp6.molarFractions[8] = 0.000746;   // н-Пентан 
    comp6.molarFractions[10] = 0.000225;   // н-Гексан
    comp6.molarFractions[11] = 0.000110; 
    comp6.molarFractions[12] = 0.000029;
    comp6.molarFractions[20] = 0.000538;





    std::cout << "Молярная масса: " << Aga8PvtCalculator::getMolarMass(comp1) << " кг/кмоль\n" << std::endl;
    
    try {
        // Расчет при p = 5 МПа, T = 250 K
        auto props = Aga8PvtCalculator::calculateFromPressure(comp1, 5, 250.0);
        props.print();
        
    } catch (const std::exception& e) {
        std::cout << "Ошибка: " << e.what() << std::endl;
    }
    
    return 0;
}