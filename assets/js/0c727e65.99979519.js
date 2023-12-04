"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[7955],{24:(e,n,i)=>{i.r(n),i.d(n,{assets:()=>o,contentTitle:()=>l,default:()=>h,frontMatter:()=>r,metadata:()=>d,toc:()=>a});var s=i(5893),t=i(1151);const r={sidebar_position:3},l="Experiment files",d={id:"user_guide/fitting/experiment_files",title:"Experiment files",description:"Description",source:"@site/docs/user_guide/fitting/experiment_files.md",sourceDirName:"user_guide/fitting",slug:"/user_guide/fitting/experiment_files",permalink:"/ChemEx/docs/user_guide/fitting/experiment_files",draft:!1,unlisted:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/user_guide/fitting/experiment_files.md",tags:[],version:"current",sidebarPosition:3,frontMatter:{sidebar_position:3},sidebar:"tutorialSidebar",previous:{title:"Starting a Fit",permalink:"/ChemEx/docs/user_guide/fitting/chemex_fit"},next:{title:"Data files",permalink:"/ChemEx/docs/user_guide/fitting/data_files"}},o={},a=[{value:"Description",id:"description",level:2},{value:"Example",id:"example",level:2},{value:"Sections",id:"sections",level:2},{value:"<code>[experiment]</code>",id:"experiment",level:3},{value:"<code>[conditions]</code>",id:"conditions",level:3},{value:"<code>label</code>",id:"label",level:4},{value:"<code>[data]</code>",id:"data",level:3},{value:"<code>error</code>",id:"error",level:4},{value:"<code>filter_offsets</code>",id:"filter_offsets",level:4},{value:"<code>[data.profiles]</code>",id:"dataprofiles",level:4},{value:"sidebar_position: 3",id:"sidebar_position-3",level:2},{value:"Overview",id:"overview",level:2},{value:"Structure of Experiment Files",id:"structure-of-experiment-files",level:2},{value:"Example of an Experiment File",id:"example-of-an-experiment-file",level:3},{value:"Detailed Sections",id:"detailed-sections",level:2},{value:"<code>[experiment]</code>",id:"experiment-1",level:3},{value:"<code>[conditions]</code>",id:"conditions-1",level:3},{value:"<code>label</code>",id:"label-1",level:4},{value:"<code>[data]</code>",id:"data-1",level:3},{value:"<code>error</code>",id:"error-1",level:4},{value:"<code>filter_offsets</code>",id:"filter_offsets-1",level:4},{value:"<code>[data.profiles]</code>",id:"dataprofiles-1",level:4}];function c(e){const n={a:"a",admonition:"admonition",code:"code",h1:"h1",h2:"h2",h3:"h3",h4:"h4",hr:"hr",li:"li",ol:"ol",p:"p",pre:"pre",strong:"strong",table:"table",tbody:"tbody",td:"td",th:"th",thead:"thead",tr:"tr",ul:"ul",...(0,t.a)(),...e.components};return(0,s.jsxs)(s.Fragment,{children:[(0,s.jsx)(n.h1,{id:"experiment-files",children:"Experiment files"}),"\n",(0,s.jsx)(n.h2,{id:"description",children:"Description"}),"\n",(0,s.jsxs)(n.p,{children:["The experiment files containing the experimental details about the datasets to\nbe fitted as well as the location of the files of the profiles. These files are\nprovided to ChemEx using the ",(0,s.jsx)(n.code,{children:"-e"})," or ",(0,s.jsx)(n.code,{children:"--experiments"})," option."]}),"\n",(0,s.jsx)(n.pre,{children:(0,s.jsx)(n.code,{className:"language-shell",children:"chemex fit -e <experiment_file1> <experiment_file2> [...]\n"})}),"\n",(0,s.jsx)(n.p,{children:"Experiment files are divided in 3 different sections:"}),"\n",(0,s.jsxs)(n.ul,{children:["\n",(0,s.jsxs)(n.li,{children:[(0,s.jsx)(n.code,{children:"[experiment]"})," includes the information about the experiment."]}),"\n",(0,s.jsxs)(n.li,{children:[(0,s.jsx)(n.code,{children:"[conditions]"})," provides the sample conditions."]}),"\n",(0,s.jsxs)(n.li,{children:[(0,s.jsx)(n.code,{children:"[data]"})," contains the details about the data location, spin system assignment\nas well as methods to estimate uncertainties or filter out some measurements."]}),"\n"]}),"\n",(0,s.jsx)(n.h2,{id:"example",children:"Example"}),"\n",(0,s.jsx)(n.p,{children:"Here is an example experiment file:"}),"\n",(0,s.jsx)(n.pre,{children:(0,s.jsx)(n.code,{className:"language-toml",metastring:'title="experiment.toml"',children:'[experiment]\nname         = "cest_15n"\ntime_t1      = 0.5\ncarrier      = 118.987\nb1_frq       = 26.3\n\n[conditions]\nh_larmor_frq = 499.243\n# p_total     = 2.0e-3\n# sample      = "A39G FF domain"\n# temperature = 1.0\n\n[data]\npath           = "../Data/26Hz/"\nerror          = "scatter"\nfilter_offsets = [[0.0, 26.0]]\n\n  [data.profiles]\n  13N = "13N-HN.out"\n  26N = "26N-HN.out"\n  28N = "28N-HN.out"\n  29N = "29N-HN.out"\n'})}),"\n",(0,s.jsx)(n.admonition,{type:"note",children:(0,s.jsxs)(n.p,{children:["The list of available key-value pairs depends on each experiment. For each\nexperiment available in ChemEx, a sample configuration file is provided in the\nsection ",(0,s.jsx)(n.a,{href:"/docs/experiments",children:"Experiments"}),". Detailed descriptions about the\nmeaning of each key are included as comments in sample config files."]})}),"\n",(0,s.jsx)(n.h2,{id:"sections",children:"Sections"}),"\n",(0,s.jsx)(n.h3,{id:"experiment",children:(0,s.jsx)(n.code,{children:"[experiment]"})}),"\n",(0,s.jsxs)(n.p,{children:["The ",(0,s.jsx)(n.code,{children:"[experiment]"})," section contains the details about the pulse sequence type\nand settings. Below are some of the keys that are commonly found in the\n",(0,s.jsx)(n.code,{children:"[experiment]"})," section:"]}),"\n",(0,s.jsxs)(n.table,{children:[(0,s.jsx)(n.thead,{children:(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.th,{children:"Name"}),(0,s.jsx)(n.th,{children:"Description"})]})}),(0,s.jsxs)(n.tbody,{children:[(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.td,{children:(0,s.jsx)(n.code,{children:"name"})}),(0,s.jsx)(n.td,{children:"Pulse sequence name."})]}),(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.td,{children:(0,s.jsx)(n.code,{children:"carrier"})}),(0,s.jsx)(n.td,{children:"Position of the carrier during the experiment, in ppm."})]}),(0,s.jsxs)(n.tr,{children:[(0,s.jsxs)(n.td,{children:[(0,s.jsx)(n.code,{children:"time_t2"}),", ",(0,s.jsx)(n.code,{children:"time_t1"})]}),(0,s.jsx)(n.td,{children:"Transverse and longitudinal relaxation delay in second (e.g. CPMG relaxation dispersion experiments)."})]}),(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.td,{children:(0,s.jsx)(n.code,{children:"pw90"})}),(0,s.jsx)(n.td,{children:"90 degree pulse width, in second."})]}),(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.td,{children:(0,s.jsx)(n.code,{children:"b1_frq"})}),(0,s.jsx)(n.td,{children:"B1 radio-frequency field strength, in Hz."})]}),(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.td,{children:(0,s.jsx)(n.code,{children:"observed_state"})}),(0,s.jsxs)(n.td,{children:["The ID of the state that is observed (one of ",(0,s.jsx)(n.code,{children:"a"}),", ",(0,s.jsx)(n.code,{children:"b"}),", ",(0,s.jsx)(n.code,{children:"c"})," or ",(0,s.jsx)(n.code,{children:"d"}),")."]})]})]})]}),"\n",(0,s.jsx)(n.h3,{id:"conditions",children:(0,s.jsx)(n.code,{children:"[conditions]"})}),"\n",(0,s.jsxs)(n.p,{children:["The section ",(0,s.jsx)(n.code,{children:"[conditions]"})," provides information about the experimental and\nsample conditions (Larmor frequency, sample temperature, protein concentration\nand labeling):"]}),"\n",(0,s.jsxs)(n.table,{children:[(0,s.jsx)(n.thead,{children:(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.th,{children:"Name"}),(0,s.jsx)(n.th,{children:"Description"})]})}),(0,s.jsxs)(n.tbody,{children:[(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.td,{children:(0,s.jsx)(n.code,{children:"h_larmor_frq"})}),(0,s.jsx)(n.td,{children:"Larmor frequency, in MHz."})]}),(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.td,{children:(0,s.jsx)(n.code,{children:"temperature"})}),(0,s.jsx)(n.td,{children:"Sample temperature, in \xb0C (optional, required by some the kinetic models)."})]}),(0,s.jsxs)(n.tr,{children:[(0,s.jsxs)(n.td,{children:[(0,s.jsx)(n.code,{children:"p_total"}),", ",(0,s.jsx)(n.code,{children:"l_total"})]}),(0,s.jsx)(n.td,{children:"Protein and ligand concentration, respectively, in M (optional, required by some the kinetic models)."})]}),(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.td,{children:(0,s.jsx)(n.code,{children:"label"})}),(0,s.jsx)(n.td,{children:"Labeling scheme of the sample."})]})]})]}),"\n",(0,s.jsxs)(n.p,{children:["While you should always provide the Larmor frequency (",(0,s.jsx)(n.code,{children:"h_larmor_frq"}),"), the\n",(0,s.jsx)(n.code,{children:"temperature"}),", ",(0,s.jsx)(n.code,{children:"p_total"})," and ",(0,s.jsx)(n.code,{children:"l_total"})," keys are only required depending on which\n",(0,s.jsx)(n.a,{href:"/ChemEx/docs/user_guide/fitting/kinetic_models",children:"kinetic model"})," is employed."]}),"\n",(0,s.jsx)(n.h4,{id:"label",children:(0,s.jsx)(n.code,{children:"label"})}),"\n",(0,s.jsxs)(n.p,{children:["The ",(0,s.jsx)(n.code,{children:"label"})," key is a list providing information about the labeling scheme of the\nsample. The labeling information is used in some experiments to take into\naccount some isotopic effects that may affect relaxation rates or the presence\nof scalar couplings."]}),"\n",(0,s.jsxs)(n.ul,{children:["\n",(0,s.jsxs)(n.li,{children:["Use ",(0,s.jsx)(n.code,{children:'"2H"'})," for deuterated samples to obtain accurate initial estimates of\nrelaxation rates based on model-free parameters."]}),"\n",(0,s.jsxs)(n.li,{children:["Use ",(0,s.jsx)(n.code,{children:'"13C"'})," for uniformly ",(0,s.jsx)("sup",{children:"13"}),"C-labeled samples to account for scalar\ncouplings in CEST experiments properly."]}),"\n"]}),"\n",(0,s.jsxs)(n.p,{children:["For example, for a uniformly ",(0,s.jsx)("sup",{children:"13"}),"C-labeled, perdeuterated sample, use\nthe following configuration:"]}),"\n",(0,s.jsx)(n.pre,{children:(0,s.jsx)(n.code,{className:"language-toml",children:'label = ["13C", "2H"]\n'})}),"\n",(0,s.jsx)(n.h3,{id:"data",children:(0,s.jsx)(n.code,{children:"[data]"})}),"\n",(0,s.jsxs)(n.p,{children:["The section ",(0,s.jsx)(n.code,{children:"[data]"})," contains the details about the data location, the\nspin-system assignments, and methods to estimate uncertainties and/or filter out\nsome data points."]}),"\n",(0,s.jsxs)(n.table,{children:[(0,s.jsx)(n.thead,{children:(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.th,{children:"Name"}),(0,s.jsx)(n.th,{children:"Description"})]})}),(0,s.jsxs)(n.tbody,{children:[(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.td,{children:(0,s.jsx)(n.code,{children:"path"})}),(0,s.jsx)(n.td,{children:"Path to the directory containing the data files."})]}),(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.td,{children:(0,s.jsx)(n.code,{children:"error"})}),(0,s.jsx)(n.td,{children:"The method for error estimation."})]}),(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.td,{children:(0,s.jsx)(n.code,{children:"filter_offsets"})}),(0,s.jsx)(n.td,{children:"The list of offsets to exclude from the calculation (in CEST experiments)."})]}),(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.td,{children:(0,s.jsx)(n.code,{children:"filter_planes"})}),(0,s.jsx)(n.td,{children:"The list of planes from which data points should be excluded from the calculation (in CEST and CPMG experiments)."})]}),(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.td,{children:(0,s.jsx)(n.code,{children:"[data.profiles]"})}),(0,s.jsx)(n.td,{children:"The subsection containing the list of the experimental profile file names are given along with their spin-system assignment"})]})]})]}),"\n",(0,s.jsx)(n.h4,{id:"error",children:(0,s.jsx)(n.code,{children:"error"})}),"\n",(0,s.jsxs)(n.p,{children:["The ",(0,s.jsx)(n.code,{children:"error"})," key is a string providing the method for error estimation."]}),"\n",(0,s.jsxs)(n.table,{children:[(0,s.jsx)(n.thead,{children:(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.th,{children:"Key value"}),(0,s.jsx)(n.th,{children:"Description"})]})}),(0,s.jsxs)(n.tbody,{children:[(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.td,{children:(0,s.jsx)(n.code,{children:'"file"'})}),(0,s.jsx)(n.td,{children:"Uncertainty values directly taken from the file, no estimation are run."})]}),(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.td,{children:(0,s.jsx)(n.code,{children:'"duplicates"'})}),(0,s.jsxs)(n.td,{children:["Estimate the uncertainty using the ",(0,s.jsx)(n.a,{href:"https://goldbook.iupac.org/html/P/P04758.html",children:"pooled standard deviation"}),". If no duplicates are found, the average of the experimental error values is returned."]})]}),(0,s.jsxs)(n.tr,{children:[(0,s.jsx)(n.td,{children:(0,s.jsx)(n.code,{children:'"scatter"'})}),(0,s.jsxs)(n.td,{children:["Estimate the uncertainty from CEST profiles. Assume the profile has additive, Gaussian noise on top of a smoothly varying signal. Adapted from ",(0,s.jsx)(n.a,{href:"https://fr.mathworks.com/matlabcentral/fileexchange/16683-estimatenoise",children:"there"}),"."]})]})]})]}),"\n",(0,s.jsx)(n.h4,{id:"filter_offsets",children:(0,s.jsx)(n.code,{children:"filter_offsets"})}),"\n",(0,s.jsx)(n.p,{children:"This is to filter out certain offsets from the CEST profiles. This can be\nuseful, for example, when decoupling sideband artifacts are present in the\nprofiles, as their locations are in general predictable. The filter is applied\nto all the profiles of the experiments."}),"\n",(0,s.jsx)(n.p,{children:"You should provide a list of pairs of values, where the first value of the pair\ncorresponds to the offset relative to the main resonance position (\u03bd) and the\nsecond value to the bandwidth around the offset (\u0394\u03bd) where points are excluded\nfrom the calculation (\u03bd \xb1 0.5 * \u0394\u03bd). All values should be given in Hz. This key\nis optional."}),"\n",(0,s.jsx)(n.pre,{children:(0,s.jsx)(n.code,{className:"language-toml",children:"filter_offsets = [\n   [0.0, 20.0],\n   [-300.0, 20.0],\n   [+300.0, 20.0],\n ]\n"})}),"\n",(0,s.jsx)(n.h4,{id:"dataprofiles",children:(0,s.jsx)(n.code,{children:"[data.profiles]"})}),"\n",(0,s.jsxs)(n.p,{children:["This subsection of the ",(0,s.jsx)(n.code,{children:"[data]"})," section contains the list of the filenames of\nthe experimental profiles and their spin-system assignments."]}),"\n",(0,s.jsxs)(n.p,{children:["The name of the spin-systems follows the Sparky-NMR peak label conventions. Each\nnucleus is designated by a ",(0,s.jsx)(n.strong,{children:"group name"}),", usually a 1- or 3-letter amino acid\nwith the position number in the protein sequence (e.g. ALA3 or A3), followed by\nan ",(0,s.jsx)(n.strong,{children:"atom name"})," (e.g. N, CA, CG1, etc.)."]}),"\n",(0,s.jsx)(n.p,{children:'For spin systems with more than one spin, use the "-" sign to separate the\ndifferent spin names. If successive spins of the system belong to the same\nresidue, the group ID (e.g. A12) can be omitted after the first spin. For\nexample, G23N-G23H and G23N-H are both valid and equivalent.'}),"\n",(0,s.jsx)(n.p,{children:"Here is an example:"}),"\n",(0,s.jsx)(n.pre,{children:(0,s.jsx)(n.code,{className:"language-toml",children:'  [data.profiles]\n  13N = "13N-HN.out"\n  26N = "26N-HN.out"\n  28N = "28N-HN.out"\n  29N = "29N-HN.out"\n'})}),"\n",(0,s.jsx)(n.admonition,{type:"note",children:(0,s.jsx)(n.p,{children:'Choose the name of the spin-systems to reflect the spin-system of interest in\neach experiment. For example, "G2N" is suitable for experiments based on a\nsingle-spin system, whereas "G2N-H" is suitable for experiments based on a\ntwo-spin system, etc. Although ChemEx has a built-in mechanism to correct\nuser-supplied spin system names automatically, we recommend that you be as\naccurate as possible when choosing them.'})}),"\n",(0,s.jsx)(n.hr,{}),"\n",(0,s.jsx)(n.h2,{id:"sidebar_position-3",children:"sidebar_position: 3"}),"\n",(0,s.jsx)(n.h1,{id:"experiment-files-1",children:"Experiment Files"}),"\n",(0,s.jsx)(n.h2,{id:"overview",children:"Overview"}),"\n",(0,s.jsxs)(n.p,{children:["Experiment files in ChemEx are used to provide detailed information about datasets and file locations for analysis. These files are accessible through the ",(0,s.jsx)(n.code,{children:"-e"})," or ",(0,s.jsx)(n.code,{children:"--experiments"})," command line option as shown below:"]}),"\n",(0,s.jsx)(n.pre,{children:(0,s.jsx)(n.code,{className:"language-shell",children:"chemex fit -e <experiment_file1> <experiment_file2> [...]\n"})}),"\n",(0,s.jsx)(n.h2,{id:"structure-of-experiment-files",children:"Structure of Experiment Files"}),"\n",(0,s.jsx)(n.p,{children:"Experiment files consist of three primary sections:"}),"\n",(0,s.jsxs)(n.ol,{children:["\n",(0,s.jsxs)(n.li,{children:[(0,s.jsx)(n.strong,{children:(0,s.jsx)(n.code,{children:"[experiment]"})}),": Contains experiment-related information."]}),"\n",(0,s.jsxs)(n.li,{children:[(0,s.jsx)(n.strong,{children:(0,s.jsx)(n.code,{children:"[conditions]"})}),": Details sample conditions."]}),"\n",(0,s.jsxs)(n.li,{children:[(0,s.jsx)(n.strong,{children:(0,s.jsx)(n.code,{children:"[data]"})}),": Includes data location, spin system assignments, and methods for uncertainty estimation or filtering measurements."]}),"\n"]}),"\n",(0,s.jsx)(n.h3,{id:"example-of-an-experiment-file",children:"Example of an Experiment File"}),"\n",(0,s.jsx)(n.pre,{children:(0,s.jsx)(n.code,{className:"language-toml",metastring:'title="experiment.toml"',children:'[experiment]\nname         = "cest_15n"\ntime_t1      = 0.5\ncarrier      = 118.987\nb1_frq       = 26.3\n\n[conditions]\nh_larmor_frq = 499.243\n# p_total     = 2.0e-3\n# sample      = "A39G FF domain"\n# temperature = 1.0\n\n[data]\npath           = "../Data/26Hz/"\nerror          = "scatter"\nfilter_offsets = [[0.0, 26.0]]\n\n  [data.profiles]\n  13N = "13N-HN.out"\n  26N = "26N-HN.out"\n  28N = "28N-HN.out"\n  29N = "29N-HN.out"\n'})}),"\n",(0,s.jsx)(n.admonition,{type:"note",children:(0,s.jsxs)(n.p,{children:["Each experiment in ChemEx comes with a sample configuration file. These files, accessible in the ",(0,s.jsx)(n.a,{href:"/docs/experiments",children:"Experiments"})," section, include detailed key-value pair descriptions."]})}),"\n",(0,s.jsx)(n.h2,{id:"detailed-sections",children:"Detailed Sections"}),"\n",(0,s.jsx)(n.h3,{id:"experiment-1",children:(0,s.jsx)(n.code,{children:"[experiment]"})}),"\n",(0,s.jsx)(n.p,{children:"This section includes pulse sequence type and settings. Common keys found here are:"}),"\n",(0,s.jsxs)(n.ul,{children:["\n",(0,s.jsxs)(n.li,{children:[(0,s.jsx)(n.code,{children:"name"}),": Name of the pulse sequence."]}),"\n",(0,s.jsxs)(n.li,{children:[(0,s.jsx)(n.code,{children:"carrier"}),": Carrier position in ppm."]}),"\n",(0,s.jsxs)(n.li,{children:[(0,s.jsx)(n.code,{children:"time_t1"}),", ",(0,s.jsx)(n.code,{children:"time_t2"}),": Relaxation delays in seconds."]}),"\n",(0,s.jsxs)(n.li,{children:[(0,s.jsx)(n.code,{children:"pw90"}),": 90 degree pulse width in seconds."]}),"\n",(0,s.jsxs)(n.li,{children:[(0,s.jsx)(n.code,{children:"b1_frq"}),": B1 field strength in Hz."]}),"\n",(0,s.jsxs)(n.li,{children:[(0,s.jsx)(n.code,{children:"observed_state"}),": Observed state ID (a, b, c, d)."]}),"\n"]}),"\n",(0,s.jsx)(n.h3,{id:"conditions-1",children:(0,s.jsx)(n.code,{children:"[conditions]"})}),"\n",(0,s.jsxs)(n.p,{children:["Provides experimental and sample conditions like Larmor frequency, sample temperature, protein concentration, and labeling. While ",(0,s.jsx)(n.code,{children:"h_larmor_frq"})," (Larmor frequency) is always required, ",(0,s.jsx)(n.code,{children:"temperature"}),", ",(0,s.jsx)(n.code,{children:"p_total"}),", and ",(0,s.jsx)(n.code,{children:"l_total"})," depend on the kinetic model used."]}),"\n",(0,s.jsx)(n.h4,{id:"label-1",children:(0,s.jsx)(n.code,{children:"label"})}),"\n",(0,s.jsxs)(n.p,{children:["The ",(0,s.jsx)(n.code,{children:"label"})," key, specifying the sample's labeling scheme, is crucial for certain experiments. Examples include ",(0,s.jsx)(n.code,{children:'"2H"'})," for deuterated samples and ",(0,s.jsx)(n.code,{children:'"13C"'})," for ",(0,s.jsx)("sup",{children:"13"}),"C-labeled samples."]}),"\n",(0,s.jsx)(n.h3,{id:"data-1",children:(0,s.jsx)(n.code,{children:"[data]"})}),"\n",(0,s.jsx)(n.p,{children:"Details data location, spin system assignments, and methods for error estimation and data filtering."}),"\n",(0,s.jsx)(n.h4,{id:"error-1",children:(0,s.jsx)(n.code,{children:"error"})}),"\n",(0,s.jsx)(n.p,{children:"Methods for error estimation include:"}),"\n",(0,s.jsxs)(n.ul,{children:["\n",(0,s.jsxs)(n.li,{children:[(0,s.jsx)(n.code,{children:'"file"'}),": Direct from file."]}),"\n",(0,s.jsxs)(n.li,{children:[(0,s.jsx)(n.code,{children:'"duplicates"'}),": Pooled standard deviation or average error values."]}),"\n",(0,s.jsxs)(n.li,{children:[(0,s.jsx)(n.code,{children:'"scatter"'}),": Gaussian noise estimation on CEST profiles."]}),"\n"]}),"\n",(0,s.jsx)(n.h4,{id:"filter_offsets-1",children:(0,s.jsx)(n.code,{children:"filter_offsets"})}),"\n",(0,s.jsx)(n.p,{children:"Used to exclude specific offsets from CEST profiles, typically provided as pairs of values indicating offset and bandwidth."}),"\n",(0,s.jsx)(n.h4,{id:"dataprofiles-1",children:(0,s.jsx)(n.code,{children:"[data.profiles]"})}),"\n",(0,s.jsx)(n.p,{children:"Lists filenames of experimental profiles and their spin-system assignments. Names should follow Sparky-NMR peak label conventions."}),"\n",(0,s.jsx)(n.admonition,{type:"note",children:(0,s.jsx)(n.p,{children:"Ensure accuracy in naming spin-systems, as it impacts experiment specificity."})})]})}function h(e={}){const{wrapper:n}={...(0,t.a)(),...e.components};return n?(0,s.jsx)(n,{...e,children:(0,s.jsx)(c,{...e})}):c(e)}},1151:(e,n,i)=>{i.d(n,{Z:()=>d,a:()=>l});var s=i(7294);const t={},r=s.createContext(t);function l(e){const n=s.useContext(r);return s.useMemo((function(){return"function"==typeof e?e(n):{...n,...e}}),[n,e])}function d(e){let n;return n=e.disableParentContext?"function"==typeof e.components?e.components(t):e.components||t:l(e.components),s.createElement(r.Provider,{value:n},e.children)}}}]);