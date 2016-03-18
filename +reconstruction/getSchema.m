function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'reconstruction', 'shan_reconstruction');
end
obj = schemaObject;
end
