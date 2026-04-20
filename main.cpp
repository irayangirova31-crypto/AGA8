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
    
    // Инициализация таблиц.
    tables::initializeTables();

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
        // Расчет при p = 5 МПа, T = 300 K
        auto props = Aga8PvtCalculator::calculateFromPressure(comp1, 5, 250.0);
        props.print();
        
    } catch (const std::exception& e) {
        std::cout << "Ошибка: " << e.what() << std::endl;
    }
    
    return 0;
}