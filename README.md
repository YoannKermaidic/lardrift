# lardrift
Code for simulating particle drift in LArTPC and computing resulting signals

# Create lardrift conda environment
`conda create --name lardriftenv`\
`conda activate lardriftenv`

# Install required libraries and dependencies
`conda install -c conda-forge pcl`           (https://pointclouds.org/) \
`conda install -c conda-forge nlohmann_json` (https://json.nlohmann.me/) \
`conda install -c conda-forge root`          (https://root.cern.ch/) \
`conda install -c conda-forge gcc`           (https://gcc.gnu.org/)


# Download and install lardrift
`cd /your/path`\
`git clone git@github.com:YoannKermaidic/lardrift.git`\
`cd /lardrift`\
`mkdir build && cd build`\
`cmake -Wno-dev ..`\
`make`\
`export PATH=/your/path/lardrift:$PATH`

# Test
`export LARDRIFT_PATH=/your/path/lardrift`\
`lardrift -cf $LARDRIFT_PATH/config/geometry.json -ef /path/to/your_COMSOL_field.txt -setup coldbox` \
For all commands: `lardrift --help`
