const mongoose = require('mongoose');

// MongoDB Schema and Model setup
const complexSchema = new mongoose.Schema({
    name: { type: String, required: true},
    description: String,
    userName: String,
    pdbStructures: [String],
    experiments: [{ type: mongoose.Schema.Types.ObjectId, ref: 'data_from_experiments' }],
    discardedPdb: [{ pdbCode: String, note: String, _id: false }],
    timestamp: String,
    active: Boolean,
});
  
// Create Mongoose Model
const ComplexData = mongoose.model('complexes', complexSchema);

module.exports = ComplexData;