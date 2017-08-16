#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/StubCleaner.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PatternGenerator.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PatternMatcher.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/MatrixBuilder.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitter.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PatternAnalyzer.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/MatrixTester.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PatternMerging.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/RoadMerging.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/NTupleMaker.h"

#include "boost/program_options.hpp"
#include <cstdlib>
#include <fstream>

// _____________________________________________________________________________
int main(int argc, char **argv) {

    // Log the input arguments
    std::cout << Color("lgreenb") << "Command:" << EndColor() << std::endl;
    std::cout << argv[0] << " ";
    for(int i = 1; i < argc-1; ++i) {
        std::cout << argv[i] << " ";
    }
    std::cout << argv[argc-1] << std::endl << std::endl;

    // _________________________________________________________________________
    // Register the program options

    // Create a struct that holds all input arguments
    ProgramOption option;

    // Declare a group of options that will be allowed only on command line
    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
        ("version"             , "Print version")
        ("help,h"              , "Produce help message")
        ("stubCleaning,C"      , "Clean stubs and pick one unique stub per layer")
        ("bankGeneration,B"    , "Generate associative memory pattern bank")
        ("patternRecognition,R", "Run associative memory pattern recognition")
        ("matrixBuilding,M"    , "Calculate matrix constants for PCA track fitting")
        ("trackFitting,T"      , "Perform track fitting")
        ("bankAnalysis,A"      , "Analyze associative memory pattern bank")
        ("matrixTesting,U"     , "Test matrix constants for PCA track fitting")
        ("patternMerging,P"    , "Run pattern merging")
        ("roadMerging,Q"       , "Run road merging")
        ("write,W"             , "Write full ntuple")
        ("no-color"            , "Turn off colored text")
        ("timing"              , "Show timing information")
        ;

    // Declare a group of options that will be allowed both on command line
    // and in config file
    po::options_description config("Configuration");
    config.add_options()
        ("input,i"      , po::value<std::string>(&option.input)->required(), "Specify input files")
        ("output,o"     , po::value<std::string>(&option.output)->required(), "Specify output file")
        ("bank,b"       , po::value<std::string>(&option.bankfile), "Specify pattern bank file")
        ("matrix,m"     , po::value<std::string>(&option.matrixfile), "Specify matrix constants file")
        ("roads"        , po::value<std::string>(&option.roadfile), "Specify file containing the roads")
        ("tracks"       , po::value<std::string>(&option.trackfile), "Specify file containing the tracks")

        ("verbosity,v"  , po::value<int>(&option.verbose)->default_value(1), "Verbosity level (-1 = very quiet; 0 = quiet, 1 = verbose, 2+ = debug)")
        ("speedup"      , po::value<int>(&option.speedup)->default_value(0), "Speed-up level")
        ("maxEvents,n"  , po::value<long long>(&option.maxEvents)->default_value(-1), "Specfiy max number of events")

        ("nLayers"      , po::value<unsigned>(&option.nLayers)->default_value(6), "Specify # of layers")
        ("nFakers"      , po::value<unsigned>(&option.nFakers)->default_value(0), "Specify # of fake superstrips")
        ("nDCBits"      , po::value<unsigned>(&option.nDCBits)->default_value(0), "Specify # of DC bits")

        // Trigger tower selection
        ("tower,t"      , po::value<unsigned>(&option.tower)->default_value(25), "Specify the trigger tower")

        // Superstrip definition
        ("superstrip,s" , po::value<std::string>(&option.superstrip)->default_value("f_sf:1_nz:8"),
        		"Specify the superstrip definition [type_factor1:XX_nz:YY(_nr_ZZ)] "
        		"type: fw(fixed width), proj (projctive), f (fountain), fop (f. optimised), ff (fountain-flower); factor1: ss(for fw),nx(for proj),nz (other); nr is only for flower. (default: sf1_nz8)")

        // Track fitting algorithm
        ("algo,f"       , po::value<std::string>(&option.algo)->default_value("LTF"), "Select track fitter -- PCA4: PCA fitter 4 params; PCA5: PCA fitter 5 params; ATF4: ATF fitter 4 params; ATF5: ATF fitter 5 params; LTF: Linearized track fitter (default: LTF)")

		//Hybrid, using flower ss or not. To use this, you must use -s ffX_nzX
		("flower_charge", po::value<int>(&option.flower_charge)->default_value(0), "Select charge sign (+1,-1,0(for passing all)...(USE IN COMBINATION with fountain-flower SS: ffX_nzY)")

        // MC truth
        ("minPt"        , po::value<float>(&option.minPt)->default_value(     3.0), "Specify min pt")
        ("maxPt"        , po::value<float>(&option.maxPt)->default_value(999999.0), "Specify max pt")
        ("minInvPt"     , po::value<float>(&option.minInvPt)->default_value(-999999.0), "Specify min signed 1/pt")
        ("maxInvPt"     , po::value<float>(&option.maxInvPt)->default_value( 999999.0), "Specify max signed 1/pt")
        ("minEta"       , po::value<float>(&option.minEta)->default_value(-2.5), "Specify min eta (signed)")
        ("maxEta"       , po::value<float>(&option.maxEta)->default_value( 2.5), "Specify max eta (signed)")
        ("minPhi"       , po::value<float>(&option.minPhi)->default_value(-M_PI), "Specify min phi (from -pi to pi)")
        ("maxPhi"       , po::value<float>(&option.maxPhi)->default_value( M_PI), "Specify max phi (from -pi to pi)")
        ("minVz"        , po::value<float>(&option.minVz)->default_value(-300.), "Specify min vertex z (cm)")
        ("maxVz"        , po::value<float>(&option.maxVz)->default_value( 300.), "Specify max vertex z (cm)")

        // Only for stub cleaning
        ("picky"        , po::value<int>(&option.picky)->default_value(1), "Specify picky level (default: 1)")

        // Only for stub cleaning and pattern recognition
        ("removeOverlap", po::bool_switch(&option.removeOverlap)->default_value(false), "remove stubs in overlap regions (default: false)")

        // Only for bank generation
        ("minFrequency" , po::value<int>(&option.minFrequency)->default_value(1), "Specify min frequency of a pattern to be stored or read")

        // Only for pattern matching
        ("maxPatterns"  , po::value<long int>(&option.maxPatterns)->default_value(999999999), "Specfiy max number of patterns")
        ("maxMisses"    , po::value<int>(&option.maxMisses)->default_value(1), "Specify max number of allowed misses (default: 1)")
        ("maxStubs"     , po::value<int>(&option.maxStubs)->default_value(4), "Specfiy max number of stubs per superstrip (default: 4)")
        ("maxRoads"     , po::value<int>(&option.maxRoads)->default_value(999999999), "Specfiy max number of roads per event")

        // Only for matrix building
        ("view"         , po::value<std::string>(&option.view)->default_value("XYZ"), "Specify fit view (e.g. XYZ, XY, RZ)")
        ("hitBits"      , po::value<unsigned>(&option.hitBits)->default_value(0), "Specify hit bits (0: all hit, 1: miss layer 1, ..., 6: miss layer 6)")

        // Only for track fitting
        ("maxChi2"      , po::value<float>(&option.maxChi2)->default_value(14.6), "Specify maximum reduced chi-squared")  // Marco De Mattia recommends 14.6 as a reasonable cut value
        ("CutPrincipals", po::bool_switch(&option.CutPrincipals)->default_value(false), "replace chi**2 cut with principals cut")
        ("minNdof"      , po::value<int>(&option.minNdof)->default_value(1), "Specify minimum degree of freedom")
        ("maxCombs"     , po::value<int>(&option.maxCombs)->default_value(999999999), "Specfiy max number of combinations per road")
        ("maxTracks"    , po::value<int>(&option.maxTracks)->default_value(999999999), "Specfiy max number of tracks per event")

        // Only for Duplicate Flag
        ("rmDuplicate", po::value<int>(&option.rmDuplicate)->default_value(2), "Duplicate removal option. The argument is the number of max stubs allowed to be shared between AM tracks")

        // Only for parameter-based duplicate removal
        ("rmParDuplicate", po::bool_switch(&option.rmParDuplicate)->default_value(false), "Parameter-based duplicate removal switch")

        //Only for alternative combination builder configuration
        ("oldCB", po::bool_switch(&option.oldCB)->default_value(true), "Use the old combination builder")
        ("FiveOfSix", po::bool_switch(&option.FiveOfSix)->default_value(true), "Do all 5/6 permutations of 6/6 roads in addition")
        ("PDDS", po::bool_switch(&option.PDDS)->default_value(false), "Switch on pairwise Delta Delta S combination cleaning")

        // Only for NTupleMaker
        ("no-trim"      , po::bool_switch(&option.no_trim)->default_value(false), "Do not trim ntuple branches")
        ;

    // Hidden options, will be allowed both on command line and in config file,
    // but will not be shown to the user.
    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("cfg"          , po::value<std::vector<std::string> >(), "Specify config file")
        ;

    // Further group the options
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config).add(hidden);

    po::options_description cfgfile_options;
    cfgfile_options.add(config).add(hidden);

    po::options_description visible("Allowed options");
    visible.add(generic).add(config);

    // Allow a positional option, which does not specify any name at all
    po::positional_options_description p;
    p.add("cfg", -1);


    // _________________________________________________________________________
    // Parse program options

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
        if (vm.count("help") || argc==1) {
            std::cout << visible << std::endl;
            return EXIT_SUCCESS;
        }

        if (vm.count("cfg")) {
            const std::vector<std::string>& cfgfiles = vm["cfg"].as<std::vector<std::string> >();
            for (unsigned i=0; i<cfgfiles.size(); ++i) {
                std::ifstream ifs(cfgfiles.at(i).c_str());
                po::store(po::parse_config_file(ifs, cfgfile_options), vm);
            }
        }

        po::notify(vm);
    } catch (const boost::program_options::error& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    // Exactly one of these options must be selected
    int vmcount = vm.count("stubCleaning")       +
                  vm.count("bankGeneration")     +
                  vm.count("patternRecognition") +
                  vm.count("matrixBuilding")     +
                  vm.count("trackFitting")       +
                  vm.count("bankAnalysis")       +
                  vm.count("matrixTesting")      +
                  vm.count("patternMerging")     +
                  vm.count("roadMerging")        +
                  vm.count("write")              ;
    if (vmcount != 1) {
        std::cerr << "ERROR: Must select exactly one of '-C', '-B', '-R', '-M', '-T', '-A', '-U', '-P', '-Q' or '-W'" << std::endl;
        //std::cout << visible << std::endl;
        return EXIT_FAILURE;
    }

    if (vm.count("no-color")) {
        slhcl1tt::NoColor();
    }

    if (vm.count("timing")) {
        slhcl1tt::ShowTiming();
    }

    // Update options
    if (option.maxEvents < 0)
        option.maxEvents = std::numeric_limits<long long>::max();

    option.nLayers = std::min(std::max(3u, option.nLayers), 8u);
    option.nFakers = std::min(std::max(0u, option.nFakers), 3u);
    option.nDCBits = std::min(std::max(0u, option.nDCBits), 4u);

    // Add options
    option.datadir = std::getenv("CMSSW_BASE");
    option.datadir += "/src/SLHCL1TrackTriggerSimulations/AMSimulation/data/";

    if (option.verbose>1)
        std::cout << option << std::endl << std::endl;


    // _________________________________________________________________________
    // Call the producers

    if (vm.count("stubCleaning")) {
        std::cout << Color("magenta") << "Start stub cleaning..." << EndColor() << std::endl;

        StubCleaner cleaner(option);
        int exitcode = cleaner.run();
        if (exitcode) {
            std::cerr << "An error occurred during stub cleaning. Exiting." << std::endl;
            return exitcode;
        }
        std::cout << "Stub cleaning " << Color("lgreenb") << "DONE" << EndColor() << "." << std::endl;

    } else if (vm.count("bankGeneration")) {
        std::cout << Color("magenta") << "Start pattern bank generation..." << EndColor() << std::endl;

        PatternGenerator generator(option);
        int exitcode = generator.run();
        if (exitcode) {
            std::cerr << "An error occurred during pattern bank generation. Exiting." << std::endl;
            return exitcode;
        }
        std::cout << "Pattern bank generation " << Color("lgreenb") << "DONE" << EndColor() << "." << std::endl;

    } else if (vm.count("patternRecognition")) {
        std::cout << Color("magenta") << "Start pattern recognition..." << EndColor() << std::endl;

        PatternMatcher matcher(option);
        int exitcode = matcher.run();
        if (exitcode) {
            std::cerr << "An error occurred during pattern recognition. Exiting." << std::endl;
            return exitcode;
        }
        std::cout << "Pattern recognition " << Color("lgreenb") << "DONE" << EndColor() << "." << std::endl;

    } else if (vm.count("matrixBuilding")) {
        std::cout << Color("magenta") << "Start PCA matrix building..." << EndColor() << std::endl;

        MatrixBuilder builder(option);
        int exitcode = builder.run();
        if (exitcode) {
            std::cerr << "An error occurred during PCA matrix building. Exiting." << std::endl;
            return exitcode;
        }
        std::cout << "Matrix building " << Color("lgreenb") << "DONE" << EndColor() << "." << std::endl;

    } else if (vm.count("trackFitting")) {
        std::cout << Color("magenta") << "Start track fitting..." << EndColor() << std::endl;

        TrackFitter fitter(option);
        int exitcode = fitter.run();
        if (exitcode) {
            std::cerr << "An error occurred during track fitting. Exiting." << std::endl;
            return exitcode;
        }
        std::cout << "Track fitting " << Color("lgreenb") << "DONE" << EndColor() << "." << std::endl;

    } else if (vm.count("bankAnalysis")) {
        std::cout << Color("magenta") << "Start pattern bank analysis..." << EndColor() << std::endl;

        PatternAnalyzer analyzer(option);
        int exitcode = analyzer.run();
        if (exitcode) {
            std::cerr << "An error occurred during pattern bank analysis. Exiting." << std::endl;
            return exitcode;
        }
        std::cout << "Pattern bank analysis " << Color("lgreenb") << "DONE" << EndColor() << "." << std::endl;

    } else if (vm.count("matrixTesting")) {
        std::cout << Color("magenta") << "Start PCA matrix testing..." << EndColor() << std::endl;

        MatrixTester tester(option);
        int exitcode = tester.run();
        if (exitcode) {
            std::cerr << "An error occurred during PCA matrix testing. Exiting." << std::endl;
            return exitcode;
        }
        std::cout << "PCA matrix testing " << Color("lgreenb") << "DONE" << EndColor() << "." << std::endl;

    } else if (vm.count("patternMerging")) {
        std::cout << Color("magenta") << "Start pattern merging..." << EndColor() << std::endl;

        PatternMerging merger(option);
        int exitcode = merger.run();
        if (exitcode) {
            std::cerr << "An error occurred during pattern merging. Exiting." << std::endl;
            return exitcode;
        }
        std::cout << "Pattern merging " << Color("lgreenb") << "DONE" << EndColor() << "." << std::endl;

    } else if (vm.count("roadMerging")) {
        std::cout << Color("magenta") << "Start road merging..." << EndColor() << std::endl;

        RoadMerging merger(option);
        int exitcode = merger.run();
        if (exitcode) {
            std::cerr << "An error occurred during road merging. Exiting." << std::endl;
            return exitcode;
        }
        std::cout << "Road merging " << Color("lgreenb") << "DONE" << EndColor() << "." << std::endl;

    } else if (vm.count("write")) {
        std::cout << Color("magenta") << "Start writing full ntuple..." << EndColor() << std::endl;

        NTupleMaker writer(option);
        int exitcode = writer.run();
        if (exitcode) {
            std::cerr << "An error occurred during writing full ntuple. Exiting." << std::endl;
            return exitcode;
        }
        std::cout << "Writing full ntuple " << Color("lgreenb") << "DONE" << EndColor() << "." << std::endl;

    }

    return EXIT_SUCCESS;
}
