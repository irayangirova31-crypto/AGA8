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
    comp.molarFractions[3] = 0.9500;   // CH4 - 95%
    comp.molarFractions[4] = 0.0300;   // C2H6 - 3%
    comp.molarFractions[5] = 0.0100;   // C3H8 - 1%
    comp.molarFractions[2] = 0.0080;   // N2 - 0.8%
    comp.molarFractions[1] = 0.0020;   // CO2 - 0.2%
    
    std::cout << "Состав: 95% CH4, 3% C2H6, 1% C3H8, 0.8% N2, 0.2% CO2" << std::endl;
    std::cout << "Молярная масса: " << Aga8PvtCalculator::getMolarMass(comp) << " кг/кмоль\n" << std::endl;
    
    try {
        // Расчет при p = 5 МПа, T = 300 K
        auto props = Aga8PvtCalculator::calculateFromPressure(comp, 5.0, 300.0);
        props.print();
        
    } catch (const std::exception& e) {
        std::cout << "Ошибка: " << e.what() << std::endl;
    }
    
    return 0;
}