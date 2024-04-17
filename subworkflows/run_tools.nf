include { BOOSTDIFF } from './tools/boostdiff'
include { GRNBOOST2 } from './tools/grnboost2'
include { Z_SCORE } from './tools/z_score'
include { DIFFCOEX } from './tools/diffcoex'
include { CHNET } from './tools/chnet'

workflow RUN_TOOLS {
    take:
        data
        tools

    main:
        runs = (1..params.n_runs)
        networks = []
        tools.each { tool -> 
            switch(tool) {
                case "boostdiff":
                    networks.add(BOOSTDIFF(data, runs))
                    break

                case "grnboost2":
                    networks.add(GRNBOOST2(data, runs))
                    break

                case "zscores":
                    networks.add(Z_SCORE(data))
                    break

                case "diffcoex":
                    networks.add(DIFFCOEX(data))
                    break

                case 'chnet':
                    networks.add(CHNET(data))
                    break
            }
        }
        for (int i = 0; i <= networks.size()-1; i++) {
            if (i == 0) {
                out = networks[i]
            } else {
                out = out.concat(networks[i])
            }
        }
        out = out.groupTuple()
        out.view()

    emit:
        out
}