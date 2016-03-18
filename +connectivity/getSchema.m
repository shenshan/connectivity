function obj = getSchema
persistent schemaObject

if isempty(schemaObject)
    common.getSchema;
    psy.getSchema;
    schemaObject = dj.Schema(dj.conn, 'connectivity', 'shan_connectivity');
end

obj = schemaObject;
end