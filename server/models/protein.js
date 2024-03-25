const mongoose = require('mongoose');

// MongoDB Schema and Model setup
const proteinSchema = new mongoose.Schema({
    name: { type: String, required: true},
    description: String,
    timestamp: String,
    active: Boolean,
    userId: { type: mongoose.Schema.Types.ObjectId, ref: 'users' },
    pdbStructures: [String],
    experiments: [{ type: mongoose.Schema.Types.ObjectId, ref: 'experiments' }],
    discardedPdb: [{ pdbCode: String, note: String, _id: false }],

});
  
// Create Mongoose Model
const ProteinData = mongoose.model('proteins', proteinSchema);

module.exports = ProteinData;