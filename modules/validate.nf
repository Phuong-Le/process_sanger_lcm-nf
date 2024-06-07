// adapted from script by Harry Hung
// checking input parameters


validParams = [
    with_match_normal: 'boolean',
    sample_paths: 'path',
    concordance_threshold: 'number',
    contamination_threshold_samples: 'number',
    contamination_threshold_match: 'number',
    marker_bed: 'path',
    marker_txt: 'path',
    vcfilter_config: '',
    mut_type: 'str',
    reference_genome: 'path',
    high_depth_region: 'path',
    phylogenetics: 'boolean',
    outdir: 'path',
    help: '',
    max_memory: '',
    max_cpus: '',
    max_time: ''
]

void validate(Map params) {
    // To save invalid parameters in this list
    invalidParams = []
    // To save invalid parameter values as "parameter : [value, issue]" in this map
    invalidValues = [:]

    params.each {
        key, value ->

        // If parameter is invalid, add it to invalidParams list and skip the following checks
        if (!validParams.keySet().contains(key)) {
            invalidParams.add(key)
            return
        }

        // Based on the value type of the parameter, perform the appropriate check
        switch (validParams[key]) {
            case 'int':
                if (value !instanceof Integer) {
                    invalidValues[key] = [value, 'integer value']
                }
                break
            
            case 'number':
                if (value !instanceof Number) {
                    invalidValues[key] = [value, 'numeric value']
                }
                break
            
            case 'str':
                if (value !instanceof String) {
                    invalidValues[key] = [value, 'string value']
                }
                break

            case 'path':
                File dir = new File(value)
                if (!(dir.exists() || dir.mkdirs())) {
                    invalidValues[key] = [value, 'directory path (invalid path or insufficient permissions)']
                }
                break
            
            case 'boolean':
                if (value !instanceof Boolean) {
                    invalidValues[key] = [value, 'boolean value']
                }
                break
            
            case '':
                break

            // Should only reach this statement if a new value type is added to validParams without adding its case above
            default:
                log.error("""
                    |Unknown value type \"${validParams[key]}\"
                    """.stripMargin())
                System.exit(1)
        }
    }

    // If invalidParams list or invalidValues map is not empty, log error messages and terminate the pipeline
    if (invalidParams || invalidValues) {
        log.error('The pipeline will now be terminated due to the following critical error(s):')

        if (invalidParams) {
            log.error("The following invalid option(s) were provided: --${invalidParams.join(', --')}.")
        }

        if (invalidValues) {
            invalidValues.each {
                key, values ->
                log.error("The provided value \"${values[0]}\" for option --${key} is not a valid ${values[1]}.")
            }
        }

        System.exit(1)
    }
}

// mut_type can only be "snp" or "indel"
if ( ! ["snp", "indel"].contains(params.mut_type) ) {
    log.error("mut_type can only be either snp or indel")
    System.exit(2)
}