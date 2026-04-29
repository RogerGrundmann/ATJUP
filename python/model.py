#!/usr/bin/env python

from pyatjup import Jupiter

class Model(object):
    """
    ATJUP Jupiter's cell structure in the atmosphere
    """
    def __init__(self, cfg_xml):
        self.jup = Jupiter()
        self.config_xml = cfg_xml

    def print_config_jup(self, config_xml):
        print("\n\n\n   Jupiter's cell structure in the atmosphere") 
        print("\n\n\n   jupiter atmosphere configuration file name = ", self.config_xml)
        print("   output path is           ", self.jup.output_path.decode('utf-8'))

    def run_Model_jup(self):
        print("\n   run_Model for the Jupiter-Atmosphere code prepared")
        self.jup.run()
        print("\n    successfully terminated Jupiter-Atmosphere code")
        print("\n")

jup = Model("config_atjup.xml")
jup.print_config_jup("config_atjup.xml")
jup.run_Model_jup()
