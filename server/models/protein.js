const mongoose = require('mongoose');

// MongoDB Schema and Model setup
const proteinSchema = new mongoose.Schema({
    name: { type: String, required: true},
    description: String,
    userName: String,
    pdbStructures: [String],
    experiments: [{ type: mongoose.Schema.Types.ObjectId, ref: 'experiments' }],
    discardedPdb: [{ pdbCode: String, note: String, _id: false }],
    timestamp: String,
    active: Boolean,
});
  
// Create Mongoose Model
const ProteinData = mongoose.model('proteins', proteinSchema);

module.exports = ProteinData;