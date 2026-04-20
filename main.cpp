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
    
    // Инициализация таблиц
    tables::initializeTables();
    
    // Состав природного газа
    GasComposition comp;
    comp.molarFractions[1] = 0.003;   // Азот
    comp.molarFractions[2] = 0.006;   // Диоксид углерода
    comp.molarFractions[3] = 0.965;   // Метан
    comp.molarFractions[4] = 0.018;   // Этан
    comp.molarFractions[5] = 0.0045;   // Пропан
    comp.molarFractions[6] = 0.001;   // н-Бутан
    comp.molarFractions[7] = 0.001;   // Изопутан
    comp.molarFractions[8] = 0.0003;   // н-Пентан
    comp.molarFractions[9] = 0.0005;   // Изопентан
    comp.molarFractions[10] = 0.0007;   // н-Гексан




    
    //std::cout << "Состав: 95% CH4, 3% C2H6, 1% C3H8, 0.8% N2, 0.2% CO2" << std::endl;
    std::cout << "Молярная масса: " << Aga8PvtCalculator::getMolarMass(comp) << " кг/кмоль\n" << std::endl;
    
    try {
        // Расчет при p = 5 МПа, T = 300 K
        auto props = Aga8PvtCalculator::calculateFromPressure(comp, 5.0, 250.0);
        props.print();
        
    } catch (const std::exception& e) {
        std::cout << "Ошибка: " << e.what() << std::endl;
    }
    
    return 0;
}