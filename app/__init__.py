#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from flask import Flask
from config import config
# celery
from celery import Celery

celery = Celery(__name__)
celery.config_from_object('jobs.celeryconfig')



def create_app(config_name):
    app = Flask(__name__)

    app.config.from_object(config[config_name])
    config[config_name].init_app(app)

    from .lib import lib as lib_blueprint
    app.register_blueprint(lib_blueprint,url_prefix="/shaolab/zMAP/lib")

    from .index import index as index_blueprint
    app.register_blueprint(index_blueprint,url_prefix='/shaolab/zMAP')

    from .zMap import zMap as zMap_blueprint
    app.register_blueprint(zMap_blueprint,url_prefix='/shaolab/zMAP/zMAP')

    from .reverse_zMap import reverse_zMap as reverse_zMap_blueprint
    app.register_blueprint(reverse_zMap_blueprint, url_prefix='/shaolab/zMAP/reverse-zMAP')

    from .QC import QC as QC_blueprint
    app.register_blueprint(QC_blueprint, url_prefix='/shaolab/zMAP/SampleQC')

    from .HVPannotation import HVPannotation as HVPannotation_blueprint
    app.register_blueprint(HVPannotation_blueprint,url_prefix='/shaolab/zMAP/HVPannotation')

    from .Sample_Subgroup import Sample_Subgroup as Sample_Subgroup_blueprint
    app.register_blueprint(Sample_Subgroup_blueprint,url_prefix='/shaolab/zMAP/Subgrouping')

    from .AssClinical import AssClinical as AssClinical_blueprint
    app.register_blueprint(AssClinical_blueprint,url_prefix='/shaolab/zMAP/AssClinical')

    from .Surv import Surv as Surv_blueprint
    app.register_blueprint(Surv_blueprint,url_prefix='/shaolab/zMAP/Survival')

    from .AssSurv import AssSurv as AssSurv_blueprint
    app.register_blueprint(AssSurv_blueprint,url_prefix='/shaolab/zMAP/AssSurv')

    from .AssMutation import AssMutation as AssMutation_blueprint
    app.register_blueprint(AssMutation_blueprint,url_prefix='/shaolab/zMAP/AssMutation')

    from .GSVA import GSVA as GSVA_blueprint
    app.register_blueprint(GSVA_blueprint,url_prefix='/shaolab/zMAP/GSVA')

    from .network import network as network_blueprint
    app.register_blueprint(network_blueprint,url_prefix='/shaolab/zMAP/Network')

    from .manual import manual as manual_blueprint
    app.register_blueprint(manual_blueprint,url_prefix='/shaolab/zMAP/manual')

    


    return app
