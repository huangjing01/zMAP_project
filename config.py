#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
flask config file
"""

import os
import hashlib
basedir = os.path.abspath(os.path.dirname(__file__))


class Config:
    SECRET_KEY = os.environ.get('SECRET_KEY') or hashlib.new(name='md5', string='@#####@#').hexdigest()
    UPLOAD_FOLDER = os.path.join(basedir, 'uploads')
    MAX_CONTENT_LENGTH = 1000 * 1024 * 1024
    @staticmethod
    def init_app(app):
        pass


class DevelopmentConfig(Config):
    """docstring for DevelopmentConfig"""
    DEBUG = True
    AAA = 'this is DevelopmentConfig'


class TestingConfig(Config):
    """docstring for TestingConfig"""
    TESTING = True
    AAA = 'this is TestingConfig'


class ProductionConfig(Config):
    """docstring for ProductionConfig"""
    AAA = 'this is ProductionConfig'
    pass

config = {
    'development': DevelopmentConfig,
    'testing': TestingConfig,
    'production': ProductionConfig,
    'default': DevelopmentConfig
}
