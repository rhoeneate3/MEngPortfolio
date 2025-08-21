#include <stdio.h>
#include <string.h>
#include <math.h>  

// External declarations for simulation parameters
extern double POWER;         
extern double PRESSURE;        
extern double L;               
extern double TEMPERATURE;
extern double N_TURNS; 
extern double TARGET_EFFICIENCY;

// Bound variables (to be set by GUI before running simulation)
extern double MIN_POWER, MAX_POWER;
extern double MIN_PRESSURE, MAX_PRESSURE;
extern double MIN_L, MAX_L;

void run_optimisation(const char* folder_path, char mode, char sim_type, int cycles,
                    double target_efficiency, double power, double pressure, double gap, double temp, double turns,
                    double min_power,  double max_power, 
                    double min_pressure, double max_pressure,
                    double min_l, double max_l);

#ifdef __cplusplus
extern "C" {
#endif

// Main simulation controller function
__declspec(dllexport) int run_simulation(
    const char* folder_path,
    const char* opt_type,      // Optimisation type ("Power", "Pressure", etc.)
    char sim_type,             // Simulation type ('C' or 'I')
    int cycles,
    double target_efficiency, double power, double pressure, double gap_length, double temperature, double turns,
    double min_power, double max_power,
    double min_pressure, double max_pressure,
    double min_l, double max_l
) {

        // Set the global bounds here directly:
    MIN_POWER = min_power;
    MAX_POWER = max_power;
    MIN_PRESSURE = min_pressure;
    MAX_PRESSURE = max_pressure;
    MIN_L = min_l;
    MAX_L = max_l;

    printf(">> [DLL] Bounds set by run_simulation:\n");
    printf("   Power       = [%.3f, %.3f]\n", MIN_POWER, MAX_POWER);
    printf("   Pressure    = [%.3f, %.3f]\n", MIN_PRESSURE, MAX_PRESSURE);
    printf("   Gap Length  = [%.4f, %.4f]\n", MIN_L, MAX_L);

    // Validate simulation type
    if (sim_type != 'C' && sim_type != 'I') {
        printf(">> Error: Invalid simulation type '%c' - use 'C' or 'I'\n", sim_type);
        return 1;
    }

    // If Optimisation is disabled
    if (strcmp(opt_type, "None") == 0) {
        POWER = power;
        PRESSURE = pressure;
        L = gap_length;
        TEMPERATURE = temperature;
        N_TURNS = turns;
        TARGET_EFFICIENCY = target_efficiency;

        printf(">> [DLL] Running non-optimized %s simulation.\n",
            (sim_type == 'C') ? "CCP" : "ICP");

        run_optimisation(folder_path, 'N', sim_type, cycles,
                 target_efficiency, power, pressure, gap_length, temperature, turns,
                 min_power, max_power,
                 min_pressure, max_pressure,
                 min_l, max_l);

        return 0;
    }

    // Optimisation type to mode character mapping
    const struct {
        const char* name;
        char mode;
    } type_map[] = {
        {"Power",       'W'},   // Use 'W' for Power Optimisation mode
        {"Pressure",    'P'},
        {"Gap Length",  'G'}
    };

    // Resolve Optimisation mode and perform conversion if needed
    for (int i = 0; i < sizeof(type_map) / sizeof(type_map[0]); i++) {
        if (strcmp(opt_type, type_map[i].name) == 0) {

            // Assign parameters
            POWER = power;
            PRESSURE   = pressure;
            L          = gap_length;
            TEMPERATURE = temperature;
            N_TURNS = turns;
            TARGET_EFFICIENCY = target_efficiency;


            printf(">> [DLL] Starting %s simulation (%s) | Cycles: %d\n", 
                   (sim_type == 'C') ? "CCP" : "ICP", opt_type, cycles);
            printf(">> [DLL] Parameters: P_in=%.3f W, P=%.1f Pa, L=%.4f m, T=%.1f K\n", 
                   power, pressure, gap_length, temperature);

            // Run Optimisation
            run_optimisation(folder_path, type_map[i].mode, sim_type, cycles, target_efficiency, 
                             power, pressure, gap_length, temperature, turns,
                             min_power, max_power, 
                             min_pressure, max_pressure, 
                             min_l, max_l);

            return 0;
        }
    }

    // Unknown Optimisation type
    printf(">> Error: Unknown Optimisation type: %s\n", opt_type);
    return 2;
}



#ifdef __cplusplus
}
#endif