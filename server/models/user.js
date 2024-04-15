const mongoose = require('mongoose');
const bcrypt = require('bcryptjs');

const userSchema = new mongoose.Schema({
    name:{ 
      type: String, 
      required: true
    },
    orcid: { 
      type: String, 
      unique: true 
    },
    username: { 
        type: String, 
        unique: true 
    },
    password: String,
    role: {
        required: true,
        type: String,
        enum: ['visitor', 'curator'],
        default: 'visitor'
    }
  });


// Encrypt password before saving
userSchema.pre('save', async function(next) {
  if (this.isModified('password')) {
    this.password = await bcrypt.hash(this.password, 12);
  }
  next();
});

// Method to check password
userSchema.methods.isValidPassword = async function(password) {
  return await bcrypt.compare(password, this.password);
};

const UserModel = mongoose.model('users', userSchema);

module.exports = UserModel;