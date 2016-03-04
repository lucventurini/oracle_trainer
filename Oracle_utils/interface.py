#!/usr/bin/env python

import sys,re,logging,os
import abc

class Interface:

    default_logger=logging.getLogger()
    stream_handler=logging.StreamHandler()
    default_logger.addHandler(stream_handler)

    @classmethod
    def setlogger(cls):
        '''Default logger to sys.stderr, if none is provided.'''
        return cls.default_logger


    def __init__(self, logger=None,
                 dataFolder=".",
                 name='',
                 new_name='',
                 record=None,
                 namespace=None
                 ):

        if logger==None:
            self.logger=self.setlogger() #Default logger
        else:
            self.logger=logger

        self.dataFolder = dataFolder
        assert os.path.exists(self.dataFolder)
        self.name=name
        self.new_name=new_name
        self.record = record
        #Transfer attributes from the argparse Namespace (or whatever object) to the instance.
        if namespace is not None:
            for attr in namespace.__dict__:
                setattr(self, attr, getattr(namespace,attr))

        self.failed=False

    @abc.abstractmethod
    def __call__(self):
        if self.failed==True:
            self.logger.error("Failed at step {0}".format(self.step))
            return False
        else:
            self.logger.warn("Finished step {0}".format(self.step))
            return True

    @abc.abstractmethod
    def define_step(self):
        '''This method will define the step of the procedure (e.g. SNAP).'''
        pass
    


    def inherit(self, parent_instance):
        '''Quick and dirty function to inherit the parameters from another class.'''

        for attr in parent_instance.__dict__:
            setattr(self, attr, getattr(parent_instance, attr))
