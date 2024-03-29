//
// Copyright (C) 2023 Richard Olsen.
// DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
//
// This file is part of MolTwister.
//
// MolTwister is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MolTwister is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MolTwister.  If not, see <https://www.gnu.org/licenses/>.
//

#pragma once

#define NA 6.022E23                 // Avogadros constant in mol^{-1}
#define R_i 0.00831446              // Ideal gas constant in kJ/(mol*K)
#define Epsilon_0 8.8541878128E-12f // Permittivity of vacuum in C^2 N^{-1} m^{-2} = C^2 J^{-1} m^{-1}
#define Unit_e 1.602176634E-19f     // Unit charge in C
#define Conv_T 120.272              // Divide to -> reduced unit, temp
#define Conv_P 1.0E-2               // Divide to -> reduced unit, momentum
#define Conv_Peta 1.0E2             // Divide to -> reduced unit, momentum
#define Conv_t 1.0E2                // Divide to -> reduced unit, time
#define Conv_v 1.0E-2               // Divide to -> reduced unit, velocity
#define Conv_press 16388.2461       // Divide to -> reduced unit, pressure (from atm)
#define Conv_force 1.0E13f          // Divide to -> reduced unit, force (from N/mol)

// We need to convert K = 1/(4 pi Eps0) to reduced units:
// K1 = (1 / (4 * pi * Eps0)) * NA  which has units N/mol C^{-2} m^2
// K2 = (K1 / 1.0E-20) which has units N/mol C^{-2} (AA)^2
// K3 = K2 * Unit_e*Unit_e which has units N/mol (AA)^2
// K = K2 / Conv_force which converts N/mol to reduced units. Calculating this we find K = 1389.32
#define Coulomb_K 1389.32f
