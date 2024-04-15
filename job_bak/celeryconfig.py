#!/usr/bin/env python3

# Broker and Backend
broker_url = 'redis://localhost:6379/0'
result_backend = 'redis://localhost:6379/1'

# import
imports = ('jobs.periodic','jobs.zMap','jobs.reverse_zMap','jobs.reverse_zMap_cancer','jobs.pca','jobs.cluster','jobs.enrichment','jobs.network')


task_serializer = 'msgpack'

result_serializer = 'json'

accept_content = ['json', 'msgpack']

result_expires = 7 * 24 * 60 * 60
# Timezone
timezone = 'Asia/Shanghai'
enable_utc = True


from datetime import timedelta

