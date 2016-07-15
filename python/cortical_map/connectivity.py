import datajoint as dj

schema = dj.schema('shan_connectivity', locals())

schema.spawn_missing_classes()