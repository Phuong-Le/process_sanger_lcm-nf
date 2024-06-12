// Courtesy of Harry Hung, original script can be found at https://github.com/sanger-bentley-group/gps-pipeline/blob/master/modules/singularity.nf
// Check if Singularity images are already pulled, otherwise pull non-existing images one by one
void singularityPreflight(LinkedHashMap workflowContainer, String singularityCacheDir) {
    log.info("Checking if all the Singularity images are available at ${singularityCacheDir}\n")

    // Get names of all images
    containers = workflowContainer.collect { it.value } as Set

    // Create the directory for saving images if not yet existed
    File cacheDir = new File(singularityCacheDir)
    cacheDir.exists() || cacheDir.mkdirs()

    // Get images that needs to be downloaded
    toDownload = [] as Set
    containers.each { container ->
        targetName = container.replace(':', '-').replace('/', '-') + '.img'
        targetFile = new File (singularityCacheDir + File.separator + targetName)
        if (!targetFile.exists()) {
            toDownload.add([container, targetName])
        }
    }

    // Download all the images that do not exist yet
    toDownload.each { container, targetName ->
        log.info("${container} is not found in cached directory ${singularityCacheDir}. Pulling now...")
        process = "singularity pull --dir ${singularityCacheDir} ${targetName} docker://${container}".execute()
        process.waitFor()

        if (process.exitValue()) {
            def errorMessage = new BufferedReader(new InputStreamReader(process.getErrorStream())).getText()

            log.info(
                """
                |Singularity Error Messages:
                |${errorMessage}
                | 
                |${container} cannot be pulled successfully. Resolve the above error and re-run the pipeline.
                |
                """.stripMargin()
            )
            System.exit(1)
        }

        log.info("${container} is pulled and saved as ${targetName}\n")
    }

    log.info("All images are ready. The workflow will resume.\n")
}