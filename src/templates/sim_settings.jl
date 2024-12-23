
"""
Generate Settings
"""

settings_tmp = mt"""
{
  "CORE": {
    "EXPERIMENT": "my_test",
    "RESULTS_DIR": "my_test_results",
    "TOTAL_TIME": 3600,
    "DELTA": 0.4,
    "MIN_DELTA_STEPS": 5,
    "MAXIMAL_REACTION_RADIUS": 0.1,
    "CHECKPOINT_FREQUENCY": 30
  },
  "EXPANDED_CELL_COMPLEX": {
    "CELL_NX": 3,
    "CELL_NY": 3,
    "CELL_DX": 1.66,
    "CELL_DY": 1.66,
    "GHOSTED": false
  },
  "INITIALIZATION": {
    "NUM_MT": 0,
    "MT_MIN_SEGMENT_INIT": 0.005,
    "MT_MAX_SEGMENT_INIT": 0.01
  },
  "GRAMMAR": {
    "GROWTH_RULES": {
      "ENABLE_GROWTH": true,
      "WITH_GROWTH_RATE_FACTOR": 100,
      "GROWTH_VELOCITY": 0.0615,
      "DIV_LENGTH": 0.075,
      "LENGTH_DIV_FACTOR": 1.2,
      "ENABLE_WOBBLE": false,
      "WOBBLE_ANGLE": 8
    },
    "RETRACTION_RULES": {
      "ENABLE_RETRACTION": true,
      "WITH_RETRACTION_RATE_FACTOR": 10,
      "RETRACTION_VELOCITY": 0.00883,
      "DIV_LENGTH_RETRACT": 0.0025
    },
    "BOUNDARY_RULES": {
      "ENABLE_STANDARD": true,
      "STANDARD_RATE": 40000,
      "ENABLE_CLASP": true,
      "CLASP_RATE": 40000
    },
    "CATASTROPHE_RULES": {
      "ENABLE_INTERMEDIATE_CIC": true,
      "INTERMEDIATE_CIC_RATE": 4000,
      "ENABLE_POSITIVE_CIC": true,
      "POSITIVE_CIC_RATE": 1000,
      "ENABLE_NEGATIVE_CIC": true,
      "NEGATIVE_CIC_RATE": 4000,
      "CATASTROPHE_ANGLE": 89
    },
    "ZIPPERING_RULES": {
      "ENABLE_ZIPPERING": true,
      "ZIPPERING_HIT_RATE": 4000,
      "ZIPPERING_GUARD_RATE": 4000,
      "ZIPPERING_RETRACTION_RATE": 10,
      "CRITICAL_ANGLE": 50,
      "SEPARATION_DISTANCE": 0.025
    },
    "CROSSOVER_RULES": {
      "ENABLE_CROSSOVER": true,
      "CROSSOVER_RATE": 40,
      "CROSSOVER_ANGLE": 50,
      "ENABLE_UNCROSSOVER": true,
      "UNCROSSOVER_RATE": 0.01
    },
    "CLASP_RULES": {
      "CLASP_ENABLE_ENTRY": false,
      "CLASP_ENTRY_RATE": 0.0005,
      "CLASP_ENTRY_ANGLE": 15,
      "CLASP_ENABLE_EXIT": false,
      "CLASP_EXIT_RATE": 40000,
      "CLASP_EXIT_ANGLE": 15,
      "CLASP_ENABLE_CAT": true,
      "CLASP_CAT_RATE": 1000,
      "CLASP_ENABLE_DETACHMENT": false,
      "CLASP_DETACHMENT_RATE": 0.01
    },
    "DESTRUCTION_RULES": {
      "ENABLE_MT_DESTRUCTION": true,
      "MT_DESTRUCTION_RATE": 10
    },
    "CREATION_RULES": {
      "ENABLE_CREATION": true,
      "CREATION_RATE": 0.0026,
      "CREATION_FACTOR": 1
    },
    "RECOVERY_RULES": {
      "ENABLE_RECOVERY": true,
      "RECOVERY_RATE": 0.016,
      "RECOVERY_FACTOR": 1
    }
  },
  "RHO_TEST_RATE": 10,
  "SIGMOID_K": 10,
  "COLLISION_DISTANCE": 0.025
}
"""

function generate_settings(data)
    Mustache.render(settings_tmp, data)
end
