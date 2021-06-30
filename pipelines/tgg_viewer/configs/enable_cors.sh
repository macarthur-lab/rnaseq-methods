gsutil cors set <( echo '[{ "origin": ["*"], "method": ["GET"] }]' ) gs://tgg-viewer-configs 
gsutil cors get gs://tgg-viewer-configs
